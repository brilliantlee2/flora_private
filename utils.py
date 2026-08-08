from collections import namedtuple, defaultdict, Counter
import gzip
import pandas as pd
from tqdm import tqdm
import numpy as np 
from matplotlib import pyplot as plt
import zipfile
import io
from fast_edit_distance import edit_distance, sub_edit_distance
from io import StringIO
import multiprocessing as mp
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import sys
import networkx as nx

# a light class for a read in fastq file
read_tuple = namedtuple('read_tuple', ['id', 'seq', 'q_letter'])
def fastq_parser(file_handle):
    while True:
        id = next(file_handle, None)
        if id is None:
            break
        seq = next(file_handle)
        next(file_handle) # skip  '+'
        q_letter = next(file_handle)
        yield read_tuple(id[1:].split()[0], seq.strip(), q_letter.strip()) #每次yield一条read的信息

# split any iterator in to batches  
def batch_iterator(iterator, batch_size):
    """generateor of batches of items in a iterator with batch_size.
    """
    batch = []
    i=0
    for entry in iterator:
        i += 1
        batch.append(entry)
        
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):  #保证批次处理的时候，最后一批不满足batch_size的那些数据，也可以被yield
        yield batch

def read_batch_generator(fastq_fns, batch_size):   #输出batch size read info
    """Generator of barches of reads from list of fastq files

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 100.
    """
    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq = fastq_parser(handle)
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
        else:
            with open(fn, "r") as handle:
                fastq = fastq_parser(handle)
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch

def reverse_complement(seq):
    '''
    Args: <str>
        queried seq
    Returns: <str>
        reverse_complement seq
    '''
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                    'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    letters = \
        [comp[base] if base in comp.keys() else base for base in seq]
    return ''.join(letters)[::-1]

def polyA_trimming_idx(seq, seed="AAAA", window=10, min_A=7, min_tail_len=8):
    """
    从 read 末端往前检测 polyA，返回 polyA 起始的绝对坐标（0-based）。
    若未检测到则返回 None。
    """
    s = seq.upper()
    anchor = s.rfind(seed)  # 最右侧 seed 的起点（绝对坐标）
    if anchor == -1:
        return None

    polyA_start = anchor  # 先把 polyA 起点放在 seed 起点
    i = anchor - 1        # 从 seed 之前的碱基开始，向左延伸

    while i >= 0:
        if s[i] == 'A':
            polyA_start = i
            i -= 1
            continue
        left = max(0, i - window + 1)
        if s[left:i+1].count('A') >= min_A:
            polyA_start = i
            i -= 1
            continue
        break

    # 确保 polyA 足够长
    if len(s) - polyA_start < min_tail_len:
        return None
    return polyA_start

def polyA_trimming_idx_neg(seq, **kwargs):
    idx_abs = polyA_trimming_idx(seq, **kwargs)  # 用上面的绝对坐标函数
    if idx_abs is None:
        return None
    return idx_abs - len(seq)  # 负数：从末尾往前的偏移

def find_pos(s, sub):
    return s.find(sub)         # 找不到就是 -1


def rfind_with_negative(s, sub):
    pos = s.rfind(sub)
    if pos == -1:
        return -1  # not found 
    return pos - len(s)

def default_count_threshold_calculation(count_array, exp_cells):
    top_count = np.sort(count_array)[::-1][:exp_cells]
    return np.quantile(top_count, 0.95)/20

def knee_plot(counts, threshold=None, out_fn = 'knee_plot.png'):
    """
    Plot knee plot using the high-confidence putative BC counts

    Args:
        counts (list): high-confidence putative BC counts
        threshold (int, optional): a line to show the count threshold. Defaults to None.
    """
    counts = sorted(counts)[::-1]
    plt.figure(figsize=(8, 8))
    plt.title(f'Barcode rank plot (from high-quality putative BC)')
    plt.loglog(counts,marker = 'o', linestyle="", alpha = 1, markersize=6)
    plt.xlabel('Barcodes')
    plt.ylabel('Read counts')
    plt.axhline(y=threshold, color='r', linestyle='--', label = 'cell calling threshold')
    plt.legend()
    plt.savefig(out_fn)

def get_bc_whitelist(raw_bc_count, full_bc_whitelist=None, exp_cells=None, out_plot_fn = None,empty_max_count = np.inf, DEFAULT_EMPTY_DROP_MIN_ED=None, DEFAULT_EMPTY_DROP_NUM=None, reverse_complement_whitelist=True):
    percentile_count_thres = default_count_threshold_calculation
    whole_whitelist = []
    def _orient(seq):
        seq = seq.strip()
        return reverse_complement(seq) if reverse_complement_whitelist else seq

    if full_bc_whitelist.endswith('.zip'):
        with zipfile.ZipFile(full_bc_whitelist) as zf:
            # check if there is only 1 file
            assert len(zf.namelist()) == 1

            with io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding="utf-8") as f:
                for line in f:
                    whole_whitelist.append(_orient(line))
    elif full_bc_whitelist.endswith(".gz"):
        with gzip.open(full_bc_whitelist, "rt", encoding="utf-8") as f:
            for line in f:
                whole_whitelist.append(_orient(line))
    else:
        with open(full_bc_whitelist, 'r') as f:
            for line in f:
                whole_whitelist.append(_orient(line))
    
    whole_whitelist = set(whole_whitelist)
    raw_bc_count = {k:v for k,v in  raw_bc_count.items() if k in whole_whitelist}
    if not raw_bc_count:
        knee_plot([], 0, out_plot_fn)
        return {}, []
    #print(len(raw_bc_count))`
    t = percentile_count_thres(list(raw_bc_count.values()), exp_cells) #t是
    knee_plot(list(raw_bc_count.values()), t, out_plot_fn)
    cells_bc = {k:v for k,v in raw_bc_count.items() if v > t}
    if not cells_bc:
        return {}, []
    
    ept_bc = []
    ept_bc_max_count = min(cells_bc.values()) #空bc最大的read支持数
    ept_bc_max_count = min(ept_bc_max_count, empty_max_count)
    #print(ept_bc_max_count)

    ept_bc_candidate = [k for k,v in raw_bc_count.items() if v < ept_bc_max_count]
    #print(len(ept_bc_candidate))
    for k in ept_bc_candidate:
        if min([edit_distance(k, x, max_ed = DEFAULT_EMPTY_DROP_MIN_ED) for x in cells_bc.keys()]) >= DEFAULT_EMPTY_DROP_MIN_ED:
            ept_bc.append(k)
            #print(len(ept_bc))
        # we don't need too much BC in this list
        if len(ept_bc) >  DEFAULT_EMPTY_DROP_NUM:
            break
    return cells_bc, ept_bc

class read_fastq:
    """This class is for mimic the Bio.SeqIO fastq record. The SeqIO is not directly used because it's slow.
    """
    def __init__(self, title, sequence, qscore, quality_map = False):
        self.id = title.split()[0].strip("@")
        self.seq = sequence
        self.qscore = qscore

def _read_and_bc_batch_generator_with_idx(fastq_fns, putative_bc_csv, batch_size):
    """Generator of barches of reads from list of fastq files with the idx of the first read
    in each batch

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 1000.
    """
    read_idx = 0
    putative_bc_f = open(putative_bc_csv, 'r')
    putative_bc_header = next(putative_bc_f)

    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in fastq_parser(handle))

                batch_iter = batch_iterator(fastq, batch_size=batch_size)
                
                for batch in batch_iter:
                    batch_len = len(batch)
                    batch_bc_df = pd.read_csv(
                        StringIO(
                            putative_bc_header + \
                            ''.join([next(putative_bc_f) for x in range(batch_len)])
                        ))
                    yield batch, read_idx, batch_bc_df
                    read_idx += batch_len
        else:
            with open(fn) as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in fastq_parser(handle))
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    batch_len = len(batch)
                    batch_bc_df = pd.read_csv(
                        StringIO( putative_bc_header + \
                        ''.join([next(putative_bc_f) for x in range(batch_len)]))
                    )

                    yield batch, read_idx, batch_bc_df #[read_id_seq_qv,....,] 0,1000,2000... df
                    read_idx += batch_len
    putative_bc_f.close()

def _match_bc_row(row, whitelist, max_ed, minQ):
    
    strand = '+'
    
    if minQ and row.putative_bc_qscore < minQ: #没有起到作用？
        return ['', '', '']

    if not row.putative_bc or row.putative_bc in whitelist: #若bc为空或者bc在白名单中，直接返回
        return [row.putative_bc, row.putative_umi, strand]
    else: #若这个barcode不为空，但是不在白名单中，需要尝试纠错
        bc = row.putative_bc
    
    best_ed = max_ed #2
    bc_hit = ''

    #对于不在white list中的barcode 这个是whitelist.csv
    #以下代码暂无问题 202508281418
    for i in whitelist:
        ed, end_idx = sub_edit_distance(i, bc, best_ed)  #直接和white list中的做比较
        if ed < best_ed:
            best_ed = ed
            bc_hit = i #最佳匹配到的一个white list中的barcode
        elif ed == best_ed:
            if not bc_hit:
                bc_hit = i
            else: 
                bc_hit = 'ambiguous'
                best_ed -= 1
                if best_ed < 0:
                    return ['', row.putative_umi, strand]
    
    if bc_hit == 'ambiguous' or bc_hit == '':
        return ['', row.putative_umi, strand] #如果矫正失败，那么将bc_corrected一列变为空，putative_umi不变
    else:
        pass
            
    out_umi = row.putative_umi
    return [bc_hit, out_umi, strand]

def assign_read_batches(r_batch,
                        whitelist_3p, whitelist_5p,
                        max_ed, gz, minQ=0,
                        emit_unmatched_fastq=True):

    """
    whitelist_3p:指的是whitelist_3p.csv
    whitelist_5p:指的是whitelist_5p.csv
    """
    read_batch, start_df_idx, df = r_batch
    df = df.fillna('')

    wl3 = set(whitelist_3p)
    wl5 = set(whitelist_5p)

    out_buffer = ''
    unmatched_fastq_buffer = ''

    # 1) 对每行同时纠错 3' 和 5'，生成 5 列
    new_cols = []
    for row in df.itertuples():
        new_cols.append(_match_bc_row_dual(row, wl3, wl5, max_ed, minQ))

    df[['BC3_corrected', 'putative_umi', 'strand',
        'BC5_corrected', 'putative_umi_5p']] = new_cols #在之前的putative上新增了这五列

    # 2) 统计“成功”的 read（你可以定义为：至少一端成功）
    ok3 = (df['BC3_corrected'] != '') & (df['putative_umi'] != '') #ok3和ok5是布尔值
    ok5 = (df['BC5_corrected'] != '') & (df['putative_umi_5p'] != '')
    demul_read_count = int((ok3 | ok5).sum()) #只有一端拆到就算成功的数量，比例应该不会大于1

    # 3) 写 fastq：至少一端 (BC+UMI) 成功才输出
    for r, bc in zip(read_batch, df.itertuples()):
        try:
            assert bc.read_id == r.id
        except AssertionError:
            err_msg("Different order in putative bc file and input fastq!", printit=True)
            sys.exit()

        side3_ok = (bc.BC3_corrected != '' and bc.putative_umi != '') #putative.csv中3‘端是否成功
        side5_ok = (bc.BC5_corrected != '' and bc.putative_umi_5p != '') #putative.csv中5‘端是否成功

        if not (side3_ok or side5_ok):
            # unmatched（沿用你原来的 A比例过滤逻辑）
            if emit_unmatched_fastq:
                putative_bc = getattr(bc, "putative_bc", "")
                if putative_bc:
                    a_ratio = putative_bc.count("A") / len(putative_bc)
                    if a_ratio <= 0.5:
                        cb3 = bc.BC3_corrected if bc.BC3_corrected else "NA"
                        ub3 = bc.putative_umi if bc.putative_umi else "NA"
                        cb5 = bc.BC5_corrected if bc.BC5_corrected else "NA"
                        ub5 = bc.putative_umi_5p if bc.putative_umi_5p else "NA"
                        #header = (f"@{cb3}_{ub3}|{cb5}_{ub5}#{bc.read_id}_{getattr(bc,'strand','+')}"
                                  #f"\tCB3:Z:{cb3}\tUB3:Z:{ub3}\tCB5:Z:{cb5}\tUB5:Z:{ub5}")
                        #header = (
                        #    f"@{bc.read_id}"
                        #    f"\tCB:Z:NA\tUB:Z:NA"
                        #    f"\tCB3:Z:{cb3}\tUB3:Z:{ub3}\tCB5:Z:{cb5}\tUB5:Z:{ub5}"
                        #    f"\tSTR:Z:{getattr(bc,'strand','+')}"
                        #)
                        header = f"@{bc.read_id}\n"
                        unmatched_fastq_buffer += header
                        unmatched_fastq_buffer += r.seq + '\n+\n' + r.qscore + '\n'
            continue

        # 4) 裁剪策略：优先用 3' 的 polyA/umi_fixed；没有就退到 5' umi_fixed（如果你有该列）
        
        has3 = (getattr(bc, "BC3_corrected", "") != "")
        has5 = (getattr(bc, "BC5_corrected", "") != "")
        
        # ======= 更稳健：只有 barcode+UMI 都有，才允许裁剪该端 =======
        trim3_ok = has3 and (getattr(bc, "putative_umi", "") != "")
        trim5_ok = has5 and (getattr(bc, "putative_umi_5p", "") != "")
        
        seq = r.seq
        qscore = r.qscore
        L = len(seq)

        umi5_start = _to_int(getattr(bc, "umi_fixed_locs_5p", None))   # 5' UMI 起点
        umi5 = getattr(bc, "putative_umi_5p", "") or ""
        start_cut = 0
        if trim5_ok and umi5_start is not None:
            start_cut = umi5_start + len(umi5)
        
        start_cut = max(0, min(L, start_cut))

        polyA = _to_int(getattr(bc, "polyA_starts", None))
        umi3_loc = _to_int(getattr(bc, "umi_fixed_locs", None))
        end_cut = L
        if trim3_ok:
            # 1) 优先 polyA（如果是负数，Python 里 end_cut = L + polyA）
            if polyA is not None:
                end_cut = (L + polyA) if polyA < 0 else polyA
            # 2) 没 polyA 用 umi_fixed_locs（沿用你原经验：cut = umi_fixed_locs - 10）
            elif umi3_loc is not None:
                cut = umi3_loc - 10
                end_cut = (L + cut) if cut < 0 else cut
            else:
                end_cut = L
        end_cut = max(0, min(L, end_cut))

        if start_cut >= end_cut:
            continue  # 剪坏/剪空了就放弃（也可以写到 unmatched 里）
        
        seq = seq[start_cut:end_cut]
        qscore = qscore[start_cut:end_cut]
        
        # 可选：太短的不输出
        if len(seq) < 30:
            continue

        # 5) header：把两端都写进去；同时给一个“主 CB/UB”（优先 3'，否则用 5'）
        """
        if side3_ok:
            CB = bc.BC3_corrected
            UB = bc.putative_umi
        else:
            CB = bc.BC5_corrected
            UB = bc.putative_umi_5p
        """

        cb3 = bc.BC3_corrected if bc.BC3_corrected else "NA"
        ub3 = bc.putative_umi if bc.putative_umi else "NA"
        cb5 = bc.BC5_corrected if bc.BC5_corrected else "NA"
        ub5 = bc.putative_umi_5p if bc.putative_umi_5p else "NA"

        if side3_ok:
            cb_main = cb3
            umi_main = ub3
        else:
            cb_main = cb5
            umi_main = ub5

        umi_main = umi_main.replace(" ", "").upper()
        if any(ch not in "ACGTN" for ch in umi_main):
            continue

        #read_name = f"{bc.read_id}_{umi_main}"
        read_name = bc.read_id

        # 你想把它们“全部输出在 readid 上”
        #read_name = f"{cb3}_{ub3}|{cb5}_{ub5}#{bc.read_id}_{bc.strand}"
        
        #out_buffer += (
        #    f"@{read_name}"
        #    f"\tCB:Z:{cb_main}\tUB:Z:{umi_main}"
        #    f"\tCB3:Z:{cb3}\tUB3:Z:{ub3}\tCB5:Z:{cb5}\tUB5:Z:{ub5}"
        #    f"\tSTR:Z:{bc.strand}\n"
        #)
        out_buffer += f"@{read_name}\n"
        out_buffer += seq + "\n+\n" + qscore + "\n"

    # 6) gzip / plain
    b_out_buffer = gzip.compress(out_buffer.encode('utf-8')) if gz else out_buffer.encode('utf-8')
    if emit_unmatched_fastq:
        b_unmatched_fastq = gzip.compress(unmatched_fastq_buffer.encode('utf-8')) if gz else unmatched_fastq_buffer.encode('utf-8')
    else:
        b_unmatched_fastq = None

    return df, b_out_buffer, demul_read_count, len(read_batch), b_unmatched_fastq



def _to_int(x):
    if x in ('', None):
        return None
    try:
        return int(float(x))  # 兼容 26.0
    except Exception:
        return None
def assign_read_batches(r_batch,
                        whitelist_3p, whitelist_5p,
                        max_ed, gz, minQ=0,
                        emit_unmatched_fastq=True):

    read_batch, start_df_idx, df = r_batch
    df = df.fillna('')

    wl3 = set(whitelist_3p)
    wl5 = set(whitelist_5p)

    out_buffer = ''
    unmatched_fastq_buffer = ''

    # 1) 对每行同时纠错 3' 和 5'，生成 5 列
    new_cols = []
    for row in df.itertuples():
        new_cols.append(_match_bc_row_dual(row, wl3, wl5, max_ed, minQ))

    df[['BC3_corrected', 'putative_umi', 'strand',
        'BC5_corrected', 'putative_umi_5p']] = new_cols #在之前的putative上新增了这五列

    # 2) 统计“成功”的 read（你可以定义为：至少一端成功）
    ok3 = (df['BC3_corrected'] != '') & (df['putative_umi'] != '') #ok3和ok5是布尔值
    ok5 = (df['BC5_corrected'] != '') & (df['putative_umi_5p'] != '')
    demul_read_count = int((ok3 | ok5).sum()) #只有一端拆到就算成功的数量，比例应该不会大于1

    # 3) 写 fastq：至少一端 (BC+UMI) 成功才输出
    for r, bc in zip(read_batch, df.itertuples()):
        try:
            assert bc.read_id == r.id
        except AssertionError:
            err_msg("Different order in putative bc file and input fastq!", printit=True)
            sys.exit()

        side3_ok = (bc.BC3_corrected != '' and bc.putative_umi != '') #putative.csv中3‘端是否成功
        side5_ok = (bc.BC5_corrected != '' and bc.putative_umi_5p != '') #putative.csv中5‘端是否成功

        if not (side3_ok or side5_ok):
            # unmatched（沿用你原来的 A比例过滤逻辑）
            if emit_unmatched_fastq:
                putative_bc = getattr(bc, "putative_bc", "")
                if putative_bc:
                    a_ratio = putative_bc.count("A") / len(putative_bc)
                    if a_ratio <= 0.5:
                        cb3 = bc.BC3_corrected if bc.BC3_corrected else "NA"
                        ub3 = bc.putative_umi if bc.putative_umi else "NA"
                        cb5 = bc.BC5_corrected if bc.BC5_corrected else "NA"
                        ub5 = bc.putative_umi_5p if bc.putative_umi_5p else "NA"
                        header = (f"@{cb3}_{ub3}|{cb5}_{ub5}#{bc.read_id}_{getattr(bc,'strand','+')}"
                                  f"\tCB3:Z:{cb3}\tUB3:Z:{ub3}\tCB5:Z:{cb5}\tUB5:Z:{ub5}")
                        unmatched_fastq_buffer += header + '\n'
                        unmatched_fastq_buffer += r.seq + '\n+\n' + r.qscore + '\n'
            continue

        # 4) 裁剪策略：优先用 3' 的 polyA/umi_fixed；没有就退到 5' umi_fixed（如果你有该列）
        
        has3 = (getattr(bc, "BC3_corrected", "") != "")
        has5 = (getattr(bc, "BC5_corrected", "") != "")
        
        # ======= 更稳健：只有 barcode+UMI 都有，才允许裁剪该端 =======
        trim3_ok = has3 and (getattr(bc, "putative_umi", "") != "")
        trim5_ok = has5 and (getattr(bc, "putative_umi_5p", "") != "")
        
        seq = r.seq
        qscore = r.qscore
        L = len(seq)

        umi5_start = _to_int(getattr(bc, "umi_fixed_locs_5p", None))   # 5' UMI 起点
        umi5 = getattr(bc, "putative_umi_5p", "") or ""
        start_cut = 0
        if trim5_ok and umi5_start is not None:
            start_cut = umi5_start + len(umi5)
        
        start_cut = max(0, min(L, start_cut))

        polyA = _to_int(getattr(bc, "polyA_starts", None))
        umi3_loc = _to_int(getattr(bc, "umi_fixed_locs", None))
        end_cut = L
        if trim3_ok:
            # 1) 优先 polyA（如果是负数，Python 里 end_cut = L + polyA）
            if polyA is not None:
                end_cut = (L + polyA) if polyA < 0 else polyA
            # 2) 没 polyA 用 umi_fixed_locs（沿用你原经验：cut = umi_fixed_locs - 10）
            elif umi3_loc is not None:
                cut = umi3_loc - 10
                end_cut = (L + cut) if cut < 0 else cut
            else:
                end_cut = L
        end_cut = max(0, min(L, end_cut))

        if start_cut >= end_cut:
            continue  # 剪坏/剪空了就放弃（也可以写到 unmatched 里）
        
        seq = seq[start_cut:end_cut]
        qscore = qscore[start_cut:end_cut]
        
        # 可选：太短的不输出
        if len(seq) < 30:
            continue

        # 5) header：把两端都写进去；同时给一个“主 CB/UB”（优先 3'，否则用 5'）
        """
        if side3_ok:
            CB = bc.BC3_corrected
            UB = bc.putative_umi
        else:
            CB = bc.BC5_corrected
            UB = bc.putative_umi_5p
        """

        cb3 = bc.BC3_corrected if bc.BC3_corrected else "NA"
        ub3 = bc.putative_umi if bc.putative_umi else "NA"
        cb5 = bc.BC5_corrected if bc.BC5_corrected else "NA"
        ub5 = bc.putative_umi_5p if bc.putative_umi_5p else "NA"

        # Keep the original read_id so downstream BAM tagging can join by read_id.
        read_name = bc.read_id

        out_buffer += f"@{read_name}\n"
        out_buffer += seq + "\n+\n" + qscore + "\n"

    # 6) gzip / plain
    b_out_buffer = gzip.compress(out_buffer.encode('utf-8')) if gz else out_buffer.encode('utf-8')
    if emit_unmatched_fastq:
        b_unmatched_fastq = gzip.compress(unmatched_fastq_buffer.encode('utf-8')) if gz else unmatched_fastq_buffer.encode('utf-8')
    else:
        b_unmatched_fastq = None

    return df, b_out_buffer, demul_read_count, len(read_batch), b_unmatched_fastq


def err_msg(msg, printit = False):
    CRED = '\033[91m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def warning_msg(msg, printit = False):
    CRED = '\033[93m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def green_msg(msg, printit = False):
    CRED = '\033[92m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,
                           pbar=True, pbar_unit='Read',pbar_func=len, 
                           schduler = 'process', *arg, **kwargs):
    """multiple processing or threading, 

    Args:
        func: function to be run parallely
        iterator: input to the function in each process/thread
        n_process (int, optional): number of cores or threads. Defaults to mp.cpu_count()-1.
        pbar (bool, optional): Whether or not to output a progres bar. Defaults to True.
        pbar_unit (str, optional): Unit shown on the progress bar. Defaults to 'Read'.
        pbar_func (function, optional): Function to calculate the total length of the progress bar. Defaults to len.
        schduler (str, optional): 'process' or 'thread'. Defaults to 'process'.

    Yields:
        return type of the func: the yield the result in the order of submit
    """
    class fake_future:
        # a fake future class to be used in single processing
        def __init__(self, rst):
            self.rst = rst
        def result(self):
            return self.rst

    if schduler == 'process':
        # make sure the number of process is not larger than the number of cores
        n_process = min(n_process-1, mp.cpu_count()-1)
        if n_process > 1:
            executor = concurrent.futures.ProcessPoolExecutor(n_process)
    elif schduler == 'thread':
        if n_process > 1:
            executor = concurrent.futures.ThreadPoolExecutor(n_process)
    else:
        green_msg('Error in multiprocessing_submit: schduler should be either process or thread', printit=True)
        sys.exit(1)

    if pbar:
        _pbar = tqdm(unit=pbar_unit, desc='Processed')
        
    # run in single process/thread if n_process < 1
    if n_process <= 1:
        for it in iterator:
            yield fake_future(func(it, *arg, **kwargs))
            if pbar:
                _pbar.update(pbar_func(it))
        return

    # A dictionary which will contain the future object
    max_queue = n_process
    futures = {}
    n_job_in_queue = 0
    
    # make sure the result is yield in the order of submit.
    job_idx = 0
    job_completed = {}

    # submit the first batch of jobs
    while n_job_in_queue < max_queue:
        i = next(iterator, None)
        if i is None:
            break
        futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
        job_idx += 1
        n_job_in_queue += 1
        job_to_yield = 0
    # yield the result in the order of submit and submit new jobs
    while True:
        # will wait until as least one job finished
        # batch size as value, release the cpu as soon as one job finished
        job = next(as_completed(futures), None)

        # yield the completed job in the order of submit  
        if job is not None:
            job_completed[futures[job][1]] = job, futures[job][0]
            del futures[job]

        # 
        if job is None and i is None and len(job_completed)==0:
            break

        # check order
        while job_to_yield in job_completed.keys():
            # update pregress bar based on batch size
            if pbar:
                _pbar.update(job_completed[job_to_yield][1])
            yield job_completed[job_to_yield][0]
            del job_completed[job_to_yield]
            
            # submit new job
            i = next(iterator, None)
            if i is not None:
                futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
                job_idx += 1
                
            job_to_yield += 1


def get_3p_features(read_info, read_ids, putative_bcs,bc_fixed_locs,putative_bc_min_qs, umis, umi_fixed_locs, post_umi_flankings, polyA_starts,read_types, BC_fixed, umi_fixed):
    part_id = read_info.id
    part_seq = read_info.seq[-30:]
    part_qv = read_info.q_letter[-30:]
    read_ids.append(part_id)
    putative_bc_min_q = None
    umi = None
    umi_fixed_loc= None
    post_umi_flanking = None
    polyA_start = None
    read_type = None
    
    BC_fixed_loc = rfind_with_negative(part_seq, BC_fixed)
    bc_fixed_locs.append(BC_fixed_loc)
    if BC_fixed_loc == -16:
        putative_bc = read_info.seq[-26:]
        putative_bcs.append(putative_bc)
        putative_bc_min_q = min([ord(x) for x in part_qv[-26:]]) -33
        putative_bc_min_qs.append(putative_bc_min_q)
        #locate umi
        find_umi_seq = read_info.seq[-36:-26] #barcode再往前10bp去找固定序列
        umi_fixed_loc_re = rfind_with_negative(find_umi_seq, umi_fixed) 
        if umi_fixed_loc_re != -1:
            
            umi_fixed_loc = umi_fixed_loc_re - 26 #相对于read的位置
            umi_fixed_locs.append(umi_fixed_loc)
            umi = read_info.seq[umi_fixed_loc - 10 : umi_fixed_loc + 5]
            umis.append(umi)
            post_umi_flanking = read_info.seq[umi_fixed_loc - 10-5 : umi_fixed_loc - 10]
            post_umi_flankings.append(post_umi_flanking)
            #鉴定polyA的起始位置
            seq_polyA = read_info.seq[umi_fixed_loc - 10 -100:umi_fixed_loc - 10]
            last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
            if last_polyA_idx: #可以检测到polyA
                polyA_start =  last_polyA_idx - 10 + umi_fixed_loc #相对于整个read
                polyA_starts.append(polyA_start)
                read_type = 1
                read_types.append(read_type)
            else:
                polyA_starts.append(polyA_start)
                read_type = 2
                read_types.append(read_type)
            
            
        else: #如果找不到umi 直接往前截取15bp 
            umi = read_info.seq[-26-15:-26] #-26-15 粗暴的认为是umi的开头位置
            umis.append(umi)
            umi_fixed_locs.append(umi_fixed_loc) #固定序列location还是为NaN
            post_umi_flanking = read_info.seq[-26-15 -5  :-26-15]
            post_umi_flankings.append(post_umi_flanking)
            #鉴定polyA的起始位置
            seq_polyA = read_info.seq[-26 -100:-26] #直接从barcode开头搜索
            last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
            if last_polyA_idx: #可以检测到polyA
                polyA_start =  last_polyA_idx + (-26) #相对于整个read
                polyA_starts.append(polyA_start)
                read_type = 3
                read_types.append(read_type)
            else:
                polyA_starts.append(polyA_start)
                read_type = 4
                read_types.append(read_type)
            
    elif BC_fixed_loc < -16: #如果barcode比较靠左
        putative_bc = read_info.seq[BC_fixed_loc-10 : BC_fixed_loc+16]
        putative_bcs.append(putative_bc)
        putative_bc_min_q = min([ord(x) for x in read_info.q_letter[BC_fixed_loc-10:BC_fixed_loc+16]]) -33
        putative_bc_min_qs.append(putative_bc_min_q)
        #locate umi
        find_umi_seq = read_info.seq[BC_fixed_loc - 10 -10 :BC_fixed_loc - 10] #barcode再往前10bp去找固定序列
        umi_fixed_loc_re = rfind_with_negative(find_umi_seq, umi_fixed)
        
            
        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re + BC_fixed_loc - 10 #相对于read的位置 修改了一个bug
            umi_fixed_locs.append(umi_fixed_loc)
            umi = read_info.seq[umi_fixed_loc - 10 : umi_fixed_loc + 5]
            umis.append(umi)
            post_umi_flanking = read_info.seq[umi_fixed_loc - 10-5 : umi_fixed_loc - 10]
            post_umi_flankings.append(post_umi_flanking)
            #if read_info.id == "250F302306011_11_158_8281_196090449_14000_1_14.58":
            #    print(BC_fixed_loc)
            #    print(umi_fixed_loc_re)
            #    print(umi)
            #    print(post_umi_flanking)
            
            #鉴定polyA的起始位置
            seq_polyA = read_info.seq[umi_fixed_loc - 10 -100 :umi_fixed_loc - 10]
            last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
            if last_polyA_idx: #可以检测到
                polyA_start =  last_polyA_idx - 10 + umi_fixed_loc #相对于整个read
                polyA_starts.append(polyA_start)
                read_type = 5
                read_types.append(read_type)
            else:
                polyA_starts.append(polyA_start)
                read_type = 6
                read_types.append(read_type)
            
        else: # 如果鉴定不到umi
            umi = read_info.seq[BC_fixed_loc - 10 -15:BC_fixed_loc - 10] #BC_fixed_loc - 10 barcode的起始
            umis.append(umi)
            umi_fixed_locs.append(umi_fixed_loc) #NaN
            post_umi_flanking = read_info.seq[BC_fixed_loc - 10 -15 -5 :BC_fixed_loc - 10-15 ]
            post_umi_flankings.append(post_umi_flanking)
            
            seq_polyA = read_info.seq[BC_fixed_loc -10 -100:BC_fixed_loc -10] #直接从barcode开头搜索
            last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
            if last_polyA_idx: #可以检测到polyA
                polyA_start =  last_polyA_idx + (BC_fixed_loc -10) #相对于整个read
                polyA_starts.append(polyA_start)
                read_type = 7
                read_types.append(read_type)
            else:
                polyA_starts.append(polyA_start)
                read_type = 8
                read_types.append(read_type)
        
        
    elif BC_fixed_loc == -1:  #鉴定不到GGAAGG，所以是-1
        #putative_bcs.append("No BC")
        putative_bcs.append(part_seq[-26:]) #直接输出后26bp
        putative_bc_min_qs.append(putative_bc_min_q)
        #虽然没有找到barcode，但是可以单独找umi
        find_umi_seq = read_info.seq[:-40] 
        umi_fixed_loc_re = rfind_with_negative(find_umi_seq, umi_fixed)  #直接在后40bp找umi
        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re #相对于read的位置
            umi_fixed_locs.append(umi_fixed_loc)
            umi = read_info.seq[umi_fixed_loc - 10 : umi_fixed_loc + 5]
            umis.append(umi)
            post_umi_flanking = read_info.seq[umi_fixed_loc - 10-5 : umi_fixed_loc - 10]
            post_umi_flankings.append(post_umi_flanking)
            #鉴定poly
            seq_polyA = read_info.seq[umi_fixed_loc - 10 -100 :umi_fixed_loc - 10]
            last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
            if last_polyA_idx: #可以检测到
                polyA_start =  last_polyA_idx - 10 + umi_fixed_loc #相对于整个read
                polyA_starts.append(polyA_start)
                read_type = 9
                read_types.append(read_type)
            else:
                polyA_starts.append(polyA_start)
                read_type = 10
                read_types.append(read_type)
        else: ### umi鉴定不到
            #putative_bcs.append(part_seq[-26:]) #直接输出后26bp
            #putative_bc_min_qs.append(putative_bc_min_q)
            umi_fixed_locs.append(umi_fixed_loc)
            umis.append(umi)
            post_umi_flankings.append(post_umi_flanking)
            polyA_starts.append(polyA_start)
            read_type = 11
            read_types.append(read_type)
                
    else:#BC_fixed_loc>-16 右半段barcode不完全 #根据umi来判断barcode，即使是残缺的barcode
        #locate umi
        find_umi_seq = read_info.seq[BC_fixed_loc - 10 -10 :BC_fixed_loc - 10] #barcode再往前10bp去找固定序列
        umi_fixed_loc_re = rfind_with_negative(find_umi_seq, umi_fixed)
        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re + BC_fixed_loc - 10 #相对于read的位置 
            umi_fixed_locs.append(umi_fixed_loc)
            umi = read_info.seq[umi_fixed_loc - 10 : umi_fixed_loc + 5]
            umis.append(umi)
            post_umi_flanking = read_info.seq[umi_fixed_loc - 10-5 : umi_fixed_loc - 10]
            post_umi_flankings.append(post_umi_flanking)
            #根据umi来定位barcode 
            putative_bc_start_loc = umi_fixed_loc + 5 
            putative_bc = read_info.seq[putative_bc_start_loc:]
            putative_bcs.append(putative_bc)
            putative_bc_min_q = min([ord(x) for x in read_info.q_letter[putative_bc_start_loc:]]) -33
            putative_bc_min_qs.append(putative_bc_min_q)

            #if read_info.id == "250F302306011_13_7768_14189_229630721_8199_2_14.36":
            #    print(BC_fixed_loc)
            #    print(umi_fixed_loc_re)
            #    print(umi)
            #    print(post_umi_flanking)
            #鉴定polyA的起始位置
            seq_polyA = read_info.seq[umi_fixed_loc - 10 -100:umi_fixed_loc - 10]
            last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
            if last_polyA_idx: #可以检测到
                polyA_start =  last_polyA_idx - 10 + umi_fixed_loc #相对于整个read
                polyA_starts.append(polyA_start)
                read_type = 12
                read_types.append(read_type)
            else:
                polyA_starts.append(polyA_start)
                read_type = 13
                read_types.append(read_type)
            
        else: #也没有找到umi
            putative_bcs.append(part_seq[-26:])
            putative_bc_min_qs.append(putative_bc_min_q)
            umi_fixed_locs.append(umi_fixed_loc)
            umis.append(umi)
            post_umi_flankings.append(post_umi_flanking)
            polyA_starts.append(polyA_start)
            read_type = 14
            read_types.append(read_type)

    return read_ids, putative_bcs, bc_fixed_locs, putative_bc_min_qs, umis, umi_fixed_locs, post_umi_flankings, polyA_starts, read_types


def get_5p_features(read_info,read_ids_5p, putative_bcs_5p, bc_fixed_locs_5p, putative_bc_min_qs_5p, umis_5p, umi_fixed_locs_5p,  BC_fixed_5p, umi_fixed_5p):
    
    part_id = read_info.id
    part_seq = read_info.seq[:30]
    part_qv = read_info.q_letter[:30]
    read_ids_5p.append(part_id)
    putative_bc_min_q = None
    umi = None
    umi_fixed_loc= None
    #post_umi_flanking = None
    polyA_start = None
    read_type = None

    BC_fixed_loc = find_pos(part_seq, BC_fixed_5p)
    bc_fixed_locs_5p.append(BC_fixed_loc)

    if BC_fixed_loc == 11: #正好在第一个11 base
        putative_bc = read_info.seq[:26]
        putative_bcs_5p.append(putative_bc)
        putative_bc_min_q = min([ord(x) for x in part_qv[:26]]) -33
        putative_bc_min_qs_5p.append(putative_bc_min_q)
        #locate umi
        find_umi_seq = read_info.seq[26:36] #barcode再往后10bp去找固定序列
        umi_fixed_loc_re = find_pos(find_umi_seq, umi_fixed_5p)

        if umi_fixed_loc_re != -1: #可以找到umi
            umi_fixed_loc = umi_fixed_loc_re + 26 #相对于read的位置
            umi_fixed_locs_5p.append(umi_fixed_loc)
            umi = read_info.seq[umi_fixed_loc : umi_fixed_loc + 15] #还是用15吧，和3‘端保持一致
            umis_5p.append(umi)

            #post_umi_flanking = read_info.seq[umi_fixed_loc + 10 + 5 : umi_fixed_loc - 10]
        else: 
            #如果找不到umi 直接往后截取15bp 
            umi = read_info.seq[26:26+15] #
            umis_5p.append(umi)
            umi_fixed_locs_5p.append(umi_fixed_loc) #固定序列location还是为NaN 两种方案，1：单独对这种read这一端的5bp固定序列做LD检测，然后往后推10bp认为是固定序列；2.先粗略找到，然后和3’端进行结合，再根据LD找到最相似的
            
    elif BC_fixed_loc > 11: #固定序列偏后，说明barcode靠后
        putative_bc = read_info.seq[BC_fixed_loc-10 : BC_fixed_loc+16]
        putative_bcs_5p.append(putative_bc)
        putative_bc_min_q = min([ord(x) for x in read_info.q_letter[BC_fixed_loc-10:BC_fixed_loc+16]]) -33
        putative_bc_min_qs_5p.append(putative_bc_min_q)
        #locate umi
        find_umi_seq = read_info.seq[BC_fixed_loc + 6 + 10  :BC_fixed_loc + 6 + 10 +10 ] #barcode再往后10bp去找固定序列
        umi_fixed_loc_re = find_pos(find_umi_seq, umi_fixed_5p)

        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re + BC_fixed_loc +6 + 10 #相对于read的位置 
            umi_fixed_locs_5p.append(umi_fixed_loc)
            umi = read_info.seq[umi_fixed_loc  : umi_fixed_loc + 15]
            umis_5p.append(umi)

        else:
            # 如果鉴定不到umi
            umi = read_info.seq[BC_fixed_loc+6+10 : BC_fixed_loc+6+10+15] #
            umis_5p.append(umi)
            umi_fixed_locs_5p.append(umi_fixed_loc)
    elif BC_fixed_loc < 11: #说明barcode靠左
        putative_bc = read_info.seq[ : BC_fixed_loc+16]
        putative_bcs_5p.append(putative_bc)
        putative_bc_min_q = min([ord(x) for x in read_info.q_letter[ : BC_fixed_loc+16]]) -33
        putative_bc_min_qs_5p.append(putative_bc_min_q)
        find_umi_seq = read_info.seq[BC_fixed_loc + 6 + 10  :BC_fixed_loc + 6 + 10 +10 ] #barcode再往后10bp去找固定序列
        umi_fixed_loc_re = find_pos(find_umi_seq, umi_fixed_5p)

        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re + BC_fixed_loc +6 + 10 #相对于read的位置 
            umi_fixed_locs_5p.append(umi_fixed_loc)
            umi = read_info.seq[umi_fixed_loc  : umi_fixed_loc + 15]
            umis_5p.append(umi)
        else:
            # 如果鉴定不到umi
            umi = read_info.seq[BC_fixed_loc+6+10 : BC_fixed_loc+6+10+15] #
            umis_5p.append(umi)
            umi_fixed_locs_5p.append(umi_fixed_loc)
    else: #找不到barcode
        #putative_bcs.append("No BC")
        putative_bcs_5p.append(part_seq[:26]) #直接输出后26bp
        putative_bc_min_qs_5p.append(putative_bc_min_q)
        
        find_umi_seq = read_info.seq[:40]  #直接在40bp扫描找umi的固定序列
        umi_fixed_loc_re = find_pos(find_umi_seq, umi_fixed_5p)  #直接在后40bp找umi
        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re
            umi_fixed_locs_5p.append(umi_fixed_loc)
            umi = read_info.seq[umi_fixed_loc  : umi_fixed_loc + 15]
            umis_5p.append(umi)
            
        else: ### umi鉴定不到
            umi_fixed_locs_5p.append(umi_fixed_loc)
            umis_5p.append(umi)

    #return read_ids, putative_bcs_5p, bc_fixed_locs_5p, putative_bc_min_qs_5p, umis_5p, umi_fixed_locs_5p

def correct_bc_one_side(bc, whitelist, max_ed):
    """return corrected_bc or '' (failed); if bc is '' -> ''"""
    if not bc:  # 空字符串 / None
        return ''
    if bc in whitelist:
        return bc

    best_ed = max_ed
    bc_hit = ''

    for i in whitelist:
        ed, _ = sub_edit_distance(i, bc, best_ed)
        if ed < best_ed:
            best_ed = ed
            bc_hit = i
        elif ed == best_ed:
            if not bc_hit:
                bc_hit = i
            else:
                bc_hit = 'ambiguous'
                best_ed -= 1
                if best_ed < 0:
                    return ''
    if bc_hit in ('', 'ambiguous'):
        return ''
    return bc_hit

def _quality_below_min(q, min_q):
    if not min_q or q is None or pd.isna(q):
        return False
    q_str = str(q).strip()
    if q_str == "" or q_str.upper() in {"NA", "NAN", "NONE"}:
        return False
    try:
        return float(q_str) < min_q
    except ValueError:
        return False

def _match_bc_row_dual(row, wl3, wl5, max_ed, minQ=0):
    strand = '+'

    # 取出两端字段（你 df 里列名按你说的）
    bc3 = getattr(row, "putative_bc", "") or ""
    umi3 = getattr(row, "putative_umi", "") or ""
    bc5 = getattr(row, "putative_bc_5p", "") or ""
    umi5 = getattr(row, "putative_umi_5p", "") or ""

    q3 = getattr(row, "putative_bc_min_qs", None)
    q5 = getattr(row, "putative_bc_min_qs_5p", None)
    if _quality_below_min(q3, minQ):
        bc3 = ""
        umi3 = ""
    if _quality_below_min(q5, minQ):
        bc5 = ""
        umi5 = ""

    # 只把“在白名单”当 ok（更符合语义）
    ok3 = (bc3 in wl3)
    ok5 = (bc5 in wl5)

    # 3 ok, 5 不 ok：纠错5
    if ok3 and (not ok5):
        bc5_corr = correct_bc_one_side(bc5, wl5, max_ed)
        umi5_out = umi5 if bc5_corr != '' else ''   # 纠错失败就清空 umi5
        return [bc3, umi3, strand, bc5_corr, umi5_out]

    # 5 ok, 3 不 ok：纠错3
    if ok5 and (not ok3):
        bc3_corr = correct_bc_one_side(bc3, wl3, max_ed)
        umi3_out = umi3 if bc3_corr != '' else ''   # 纠错失败就清空 umi3
        return [bc3_corr, umi3_out, strand, bc5, umi5]

    # 两端都不 ok：两端都纠错
    bc3_corr = correct_bc_one_side(bc3, wl3, max_ed)
    bc5_corr = correct_bc_one_side(bc5, wl5, max_ed)
    umi3_out = umi3 if bc3_corr != '' else ''
    umi5_out = umi5 if bc5_corr != '' else ''
    return [bc3_corr, umi3_out, strand, bc5_corr, umi5_out]

def assign_read(
    fastq_fns=None,
    fastq_out=None,
    putative_bc_csv=None,
    whitelsit_3p=None,
    whitelsit_5p=None,
    max_ed=None,
    n_process=None,
    batchsize=None,
    minQ=0,
    write_fastq_out=True,
    write_unmatched_fastq=True,
):
    gz = bool(fastq_out and fastq_out.endswith('.gz'))
    out_dir = os.path.dirname(fastq_out) if fastq_out else os.getcwd()
    unmatched_out = os.path.join(out_dir, "unmatched_reads.fastq.gz")
        
    r_batches = \
        _read_and_bc_batch_generator_with_idx(fastq_fns, putative_bc_csv, batchsize)
    
    whitelist_3p_list = [] 
    with open(whitelsit_3p, 'r') as f:
        for line in f:
            whitelist_3p_list.append(line.split('-')[0].strip())

    whitelist_5p_list = [] 
    with open(whitelsit_5p, 'r') as f:
        for line in f:
            whitelist_5p_list.append(line.split('-')[0].strip())

    if n_process == 1:
        demul_count_tot = 0
        count_tot = 0
        df_list = []
        output_handle = open(fastq_out, 'wb') if write_fastq_out else None
        unmatched_handle = open(unmatched_out, 'wb') if write_unmatched_fastq else None
        try:
            pbar = tqdm(unit="Reads", desc='Processed')
            for r_batch in r_batches:
                df, b_fast_str, demul_count, read_count, b_unmatched_fastq = assign_read_batches(
                    r_batch,
                    whitelist_3p_list,
                    whitelist_5p_list,
                    max_ed,
                    gz,
                    minQ=minQ,
                    emit_unmatched_fastq=write_unmatched_fastq,
                )
                demul_count_tot += demul_count
                count_tot += read_count
                if output_handle is not None:
                    output_handle.write(b_fast_str)
                if unmatched_handle is not None and b_unmatched_fastq:
                    unmatched_handle.write(b_unmatched_fastq)
                df_list.append(df)
                pbar.update(read_count) #
        finally:
            if output_handle is not None:
                output_handle.close()
            if unmatched_handle is not None:
                unmatched_handle.close()
        big_df = pd.concat(df_list, ignore_index=True) if df_list else pd.DataFrame()
        if write_fastq_out:
            green_msg(f"Reads assignment completed. Demultiplexed read saved in {fastq_out}!")
        else:
            green_msg("Reads assignment completed. Demultiplexed FASTQ writing skipped.")
        if write_unmatched_fastq:
            green_msg(f"Unmatched reads saved in: {unmatched_out}!")
        else:
            green_msg("Unmatched FASTQ writing skipped.")
        
    else:
        rst_futures = multiprocessing_submit(assign_read_batches, 
                            r_batches, 
                            n_process=n_process,
                            schduler = "process",
                            pbar_func=lambda x: len(x[0]),
                            whitelist_3p = whitelist_3p_list,
                            whitelist_5p = whitelist_5p_list,
                            max_ed = max_ed,
                            gz = gz,
                            minQ = minQ,
                            emit_unmatched_fastq=write_unmatched_fastq)

        demul_count_tot = 0
        count_tot = 0
        df_list = []
        output_handle = open(fastq_out, 'wb') if write_fastq_out else None
        unmatched_handle = open(unmatched_out, 'wb') if write_unmatched_fastq else None
        try:
            for f in rst_futures:
                df, b_fast_str, demul_count, read_count, b_unmatched_fastq = f.result()
                demul_count_tot += demul_count
                count_tot += read_count
                if output_handle is not None:
                    output_handle.write(b_fast_str)

                if unmatched_handle is not None and b_unmatched_fastq:
                    unmatched_handle.write(b_unmatched_fastq)
                # 保存df
                df_list.append(df)
        finally:
            if output_handle is not None:
                output_handle.close()
            if unmatched_handle is not None:
                unmatched_handle.close()
        big_df = pd.concat(df_list, ignore_index=True)
    
        if write_fastq_out:
            green_msg(f"Reads assignment completed. Demultiplexed read saved in {fastq_out}!")
        else:
            green_msg("Reads assignment completed. Demultiplexed FASTQ writing skipped.")
        if write_unmatched_fastq:
            green_msg(f"Unmatched reads saved in: {unmatched_out}!")
        else:
            green_msg("Unmatched FASTQ writing skipped.")
    
    return demul_count_tot, count_tot,big_df

def norm_bc(x):
    if pd.isna(x): return ""
    x = str(x).strip()
    return "" if x in ("", "NA", "nan", "None") else x

def is_missing(x):
    if pd.isna(x):
        return True
    s = str(x).strip()
    return s in ("", "NA", "nan", "None")

def norm_seq(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    return "" if s in ("", "NA", "NAN", "NONE") else s

_rc_map = str.maketrans("ACGTN", "TGCAN")

def revcomp(seq):
    seq = norm_seq(seq)
    return seq.translate(_rc_map)[::-1]

def remove_middle_6(seq, middle6):
    seq = norm_seq(seq)
    if not seq:
        return ""
    if len(seq) == 26 and seq[10:16] == middle6:
        return seq[:10] + seq[16:]
    return ""

def strip_fixed_3p(seq, middle6="GCTACC"):
    return remove_middle_6(seq, middle6)

def strip_fixed_5p(seq, middle6="CCTTCC"):
    return remove_middle_6(seq, middle6)

def build_graph_from_pair_counts(pair_counts_kept, col5="BC5n", col3="BC3n", wcol="support_reads"):
    cols = [col5, col3, wcol]
    has_umis = "support_umis" in pair_counts_kept.columns
    if has_umis:
        cols.append("support_umis")
    df = pair_counts_kept[cols].copy()
    df[col5] = df[col5].map(norm_bc)
    df[col3] = df[col3].map(norm_bc)
    df = df[(df[col5] != "") & (df[col3] != "")]
    df[wcol] = df[wcol].astype(int)
    if has_umis:
        df["support_umis"] = df["support_umis"].astype(int)

    a = df[[col5, col3]].min(axis=1)
    b = df[[col5, col3]].max(axis=1)
    uv = pd.DataFrame({"a": a, "b": b, "w": df[wcol].values})
    if has_umis:
        uv["support_umis"] = df["support_umis"].values
        edge_agg = uv.groupby(["a", "b"], as_index=False).agg(
            w=("w", "sum"),
            support_umis=("support_umis", "sum"),
        )
    else:
        edge_agg = uv.groupby(["a", "b"], as_index=False)["w"].sum()

    graph = nx.Graph()
    if has_umis:
        for a, b, w, support_umis in edge_agg.itertuples(index=False):
            graph.add_edge(a, b, weight=int(w), support_umis=int(support_umis))
    else:
        for a, b, w in edge_agg.itertuples(index=False):
            graph.add_edge(a, b, weight=int(w), support_umis=0)
    return graph, edge_agg

def component_category(graph):
    n = graph.number_of_nodes()
    m = graph.number_of_edges()
    m_self = sum(1 for u, v in graph.edges() if u == v)
    m_no_self = m - m_self

    if n == 1 and m == 1 and m_self == 1:
        return "self_only"
    if n == 2 and m_no_self == 1:
        return "pair_only"
    if n == 3 and m_no_self == 3:
        return "triangle_only"
    if n == 4 and m_no_self == 6:
        return "clique4_only"
    if n == 5 and m_no_self == 10:
        return "clique5_only"
    return "other"

def collect_components_by_type(graph):
    by_type = {}
    for comp_nodes in nx.connected_components(graph):
        subgraph = graph.subgraph(comp_nodes).copy()
        category = component_category(subgraph)
        by_type.setdefault(category, []).append(subgraph)
    return by_type

def build_core_cells(
    pair_counts_kept,
    allowed_types=("self_only", "pair_only", "triangle_only", "clique4_only", "clique5_only"),
    include_other_components=True,
    max_other_component_barcodes=8,
    col5="BC5n",
    col3="BC3n",
    wcol="support_reads",
):
    graph, edge_agg = build_graph_from_pair_counts(pair_counts_kept, col5=col5, col3=col3, wcol=wcol)
    by_type = collect_components_by_type(graph)

    core_components = []
    for cell_type in allowed_types:
        core_components.extend(by_type.get(cell_type, []))
    if include_other_components:
        for subgraph in by_type.get("other", []):
            if subgraph.number_of_nodes() <= int(max_other_component_barcodes):
                core_components.append(subgraph)

    cell_records = []
    barcode2cell = {}
    for idx, subgraph in enumerate(core_components, 1):
        cell_id = f"cell_{idx:06d}"
        barcodes = sorted(subgraph.nodes())
        cell_records.append(
            {
                "cell_id": cell_id,
                "type": component_category(subgraph),
                "n_barcodes": len(barcodes),
                "barcodes": barcodes,
            }
        )
        for barcode in barcodes:
            barcode2cell[barcode] = cell_id

    return graph, edge_agg, pd.DataFrame(cell_records), barcode2cell

def compute_top1_dominance(pair_counts_kept, col5="BC5n", col3="BC3n", wcol="support_reads"):
    pc = pair_counts_kept[[col5, col3, wcol]].copy()
    pc[col5] = pc[col5].map(norm_bc)
    pc[col3] = pc[col3].map(norm_bc)
    pc = pc[(pc[col5] != "") & (pc[col3] != "")]
    pc[wcol] = pc[wcol].astype(int)

    long = pd.concat(
        [
            pc.rename(columns={col5: "barcode", col3: "partner", wcol: "w"}),
            pc.rename(columns={col3: "barcode", col5: "partner", wcol: "w"}),
        ],
        ignore_index=True,
    )
    bp = long.groupby(["barcode", "partner"], as_index=False)["w"].sum()
    bp_sorted = bp.sort_values(["barcode", "w"], ascending=[True, False])
    top1 = bp_sorted.groupby("barcode").head(1).rename(columns={"partner": "top1_partner", "w": "top1_w"})
    sum_all = bp.groupby("barcode")["w"].sum().reset_index(name="sum_all")
    dom = sum_all.merge(top1[["barcode", "top1_partner", "top1_w"]], on="barcode", how="left")
    dom["dominance"] = dom["top1_w"] / dom["sum_all"]
    return dom

def assign_reads_to_cells(
    df_reads,
    barcode2cell,
    dom_table,
    dominance_min=0.80,
    absorb_unassigned_paired=True,
    bc5_col="BC5n",
    bc3_col="BC3n",
    read_id_col="read_id",
):
    df = df_reads.copy()
    df["_b5"] = df[bc5_col].map(norm_bc)
    df["_b3"] = df[bc3_col].map(norm_bc)

    df["has5"] = df["_b5"] != ""
    df["has3"] = df["_b3"] != ""
    df["read_kind"] = df["has5"].astype(int).astype(str) + "+" + df["has3"].astype(int).astype(str)

    df["cell_A_5"] = df["_b5"].map(barcode2cell)
    df["cell_A_3"] = df["_b3"].map(barcode2cell)
    df["cell_A"] = df["cell_A_5"].combine_first(df["cell_A_3"])
    df["cell_A_conflict"] = (
        df["cell_A_5"].notna() & df["cell_A_3"].notna() & (df["cell_A_5"] != df["cell_A_3"])
    )

    dom = dom_table.set_index("barcode") if len(dom_table) else pd.DataFrame()

    def try_absorb(bc):
        if bc == "" or len(dom) == 0 or bc not in dom.index:
            return None
        row = dom.loc[bc]
        if pd.isna(row["top1_partner"]) or float(row["dominance"]) < dominance_min:
            return None
        return barcode2cell.get(row["top1_partner"], None)

    need_B = df["cell_A"].isna() & (df["has5"] ^ df["has3"])
    df.loc[need_B & df["has5"], "cell_B"] = df.loc[need_B & df["has5"], "_b5"].map(try_absorb)
    df.loc[need_B & df["has3"], "cell_B"] = df.loc[need_B & df["has3"], "_b3"].map(try_absorb)

    if absorb_unassigned_paired:
        need_B_paired = df["cell_A"].isna() & df["has5"] & df["has3"]
        df.loc[need_B_paired, "cell_B_5"] = df.loc[need_B_paired, "_b5"].map(try_absorb)
        df.loc[need_B_paired, "cell_B_3"] = df.loc[need_B_paired, "_b3"].map(try_absorb)

        both_same = (
            need_B_paired
            & df["cell_B_5"].notna()
            & df["cell_B_3"].notna()
            & (df["cell_B_5"] == df["cell_B_3"])
        )
        only_5 = need_B_paired & df["cell_B_5"].notna() & df["cell_B_3"].isna()
        only_3 = need_B_paired & df["cell_B_5"].isna() & df["cell_B_3"].notna()
        df.loc[both_same, "cell_B"] = df.loc[both_same, "cell_B_5"]
        df.loc[only_5, "cell_B"] = df.loc[only_5, "cell_B_5"]
        df.loc[only_3, "cell_B"] = df.loc[only_3, "cell_B_3"]

    df["cell_id"] = df["cell_A"].combine_first(df["cell_B"])

    stats = {
        "n_total_reads": len(df),
        "n_paired": int((df["read_kind"] == "1+1").sum()),
        "n_single_5": int((df["read_kind"] == "1+0").sum()),
        "n_single_3": int((df["read_kind"] == "0+1").sum()),
        "n_assigned_A": int(df["cell_A"].notna().sum()),
        "n_assigned_B_only": int((df["cell_A"].isna() & df["cell_B"].notna()).sum()) if "cell_B" in df else 0,
        "n_unassigned": int(df["cell_id"].isna().sum()),
        "n_conflict_paired_A": int(df["cell_A_conflict"].sum()),
        "n_paired_absorbed_B": int((df["cell_A"].isna() & df["has5"] & df["has3"] & df["cell_B"].notna()).sum())
        if absorb_unassigned_paired
        else 0,
    }
    df = df.drop(columns=["_b5", "_b3", "cell_B_5", "cell_B_3"], errors="ignore")
    return df, stats

def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def iter_fastq(handle):
    while True:
        h = handle.readline()
        if not h:
            return
        s = handle.readline()
        p = handle.readline()
        q = handle.readline()
        if not q:
            return
        yield h, s, p, q

def extract_reads_and_filter_df_by_raw(
    df,
    raw_fastq_gz,
    out_fastq_gz,
    read_id_col="read_id",
    remove_found=True,
    return_missing_ids=False,
):
    if read_id_col not in df.columns:
        raise KeyError(f"df missing column: {read_id_col}")

    df2 = df.copy()
    df2[read_id_col] = df2[read_id_col].astype(str).str.strip()
    df2 = df2[df2[read_id_col].notna() & (df2[read_id_col] != "")].copy()

    id_set = set(df2[read_id_col].tolist())
    n_target_unique = len(id_set)
    n_target_rows = len(df2)
    found_ids = set()
    found = 0
    scanned = 0

    with open_maybe_gz(raw_fastq_gz, "rt") as fin, gzip.open(out_fastq_gz, "wt") as fout:
        for h, s, p, q in iter_fastq(fin):
            scanned += 1
            rid = h[1:].strip().split()[0]
            if rid in id_set:
                fout.write(h)
                fout.write(s)
                fout.write(p)
                fout.write(q)
                found += 1
                found_ids.add(rid)
                if remove_found:
                    id_set.remove(rid)
                    if not id_set:
                        break

    df_kept = df2[df2[read_id_col].isin(found_ids)].copy()
    stats = {
        "target_rows": int(n_target_rows),
        "target_unique_ids": int(n_target_unique),
        "found_unique_ids": int(len(found_ids)),
        "written_reads": int(found),
        "scanned_reads": int(scanned),
        "dropped_rows": int(n_target_rows - len(df_kept)),
        "missing_unique_ids": int(n_target_unique - len(found_ids)),
        "out_fastq_gz": out_fastq_gz,
    }
    if return_missing_ids:
        stats["missing_ids"] = sorted(list(set(df2[read_id_col]) - found_ids))
    return df_kept, stats

def norm_umi(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    return "" if s in ("", "NA", "NAN", "NONE") else s

def a_ratio(seq):
    if not seq:
        return 0.0
    s = seq.upper()
    return s.count("A") / len(s) if len(s) > 0 else 0.0

def resolve_pair_min(pair_counts_all, pair_min=None, auto_pair_min_floor=10, auto_pair_min_quantile=0.1):
    if pair_min is not None:
        return int(pair_min), "manual"
    if len(pair_counts_all) == 0:
        return int(auto_pair_min_floor), "auto_empty"
    q = float(np.quantile(pair_counts_all["support_reads"], auto_pair_min_quantile))
    resolved = max(int(auto_pair_min_floor), int(q))
    return resolved, "auto_quantile"


def _canonical_edge_umi_key(row):
    b5 = row["_b5"]
    b3 = row["_b3"]
    u5 = row["_umi5"]
    u3 = row["_umi3"]
    if b5 == b3:
        left, right = sorted([u5, u3])
        return f"{left}|{right}"
    if b5 <= b3:
        return f"{u5}|{u3}"
    return f"{u3}|{u5}"


def build_pair_counts_with_umis(df_paired):
    if len(df_paired) == 0:
        return pd.DataFrame(columns=["BC5n", "BC3n", "support_reads", "support_umis"])

    tmp = df_paired.copy()
    tmp["edge_umi_key"] = tmp.apply(_canonical_edge_umi_key, axis=1)
    pair_read_counts = (
        tmp.groupby(["pair_u", "pair_v"])
        .size()
        .reset_index(name="support_reads")
    )
    pair_umi_counts = (
        tmp.groupby(["pair_u", "pair_v"])["edge_umi_key"]
        .nunique()
        .reset_index(name="support_umis")
    )
    pair_counts = (
        pair_read_counts.merge(pair_umi_counts, on=["pair_u", "pair_v"], how="left")
        .rename(columns={"pair_u": "BC5n", "pair_v": "BC3n"})
        .sort_values("support_reads", ascending=False)
        .reset_index(drop=True)
    )
    pair_counts["support_reads"] = pair_counts["support_reads"].astype(int)
    pair_counts["support_umis"] = pair_counts["support_umis"].astype(int)
    return pair_counts


def _top1_map_from_metric(pc, metric_col):
    if len(pc) == 0:
        return {}
    long = pd.concat(
        [
            pc.rename(columns={"BC5n": "barcode", "BC3n": "partner", metric_col: "w"}),
            pc.rename(columns={"BC3n": "barcode", "BC5n": "partner", metric_col: "w"}),
        ],
        ignore_index=True,
    )
    top1 = (
        long.sort_values(["barcode", "w"], ascending=[True, False])
        .groupby("barcode")
        .head(1)[["barcode", "w"]]
        .reset_index(drop=True)
    )
    return top1.set_index("barcode")["w"].to_dict()


def _endpoint_pass_with_umi(edge_reads, edge_umis, top_reads, top_umis, alpha_reads, alpha_umis):
    read_rel = (edge_reads / top_reads) if top_reads else 0.0
    umi_rel = (edge_umis / top_umis) if top_umis else 0.0
    if read_rel >= alpha_reads:
        return True
    if umi_rel >= alpha_umis and read_rel >= 0.5 * alpha_reads:
        return True
    if read_rel >= 0.75 * alpha_reads and umi_rel >= 0.5 * alpha_umis:
        return True
    return False


def _component_self_and_pair_edges(subgraph):
    self_edges = []
    pair_edges = []
    for u, v, data in subgraph.edges(data=True):
        record = {
            "u": u,
            "v": v,
            "support_reads": int(data.get("weight", 0)),
            "support_umis": int(data.get("support_umis", 0)),
        }
        if u == v:
            self_edges.append(record)
        else:
            pair_edges.append(record)
    return self_edges, pair_edges


def _self_metric_map(self_edges):
    out = {}
    for edge in self_edges:
        out[edge["u"]] = edge
    return out


def _pairwise_fraction(pair_edges, self_edges, metric):
    pair_sum = sum(e[metric] for e in pair_edges)
    self_sum = sum(e[metric] for e in self_edges)
    total = pair_sum + self_sum
    return (pair_sum / total) if total else 0.0


def _weakest_to_median_ratio(edges, metric):
    vals = sorted([e[metric] for e in edges if e[metric] > 0])
    if not vals:
        return 0.0
    med = float(np.median(vals))
    if med <= 0:
        return 0.0
    return float(vals[0] / med)


def _node_pair_fraction(subgraph, node, self_map, metric):
    incident = []
    for _, nbr, data in subgraph.edges(node, data=True):
        if nbr == node:
            continue
        incident.append(int(data.get("weight" if metric == "support_reads" else "support_umis", 0)))
    pair_sum = sum(incident)
    self_val = int(self_map.get(node, {}).get(metric, 0))
    denom = pair_sum + self_val
    return (pair_sum / denom) if denom else 0.0


def _nth_incident_pair_to_self(subgraph, node, self_map, metric, nth_largest):
    vals = []
    attr = "weight" if metric == "support_reads" else "support_umis"
    for _, nbr, data in subgraph.edges(node, data=True):
        if nbr == node:
            continue
        vals.append(int(data.get(attr, 0)))
    vals = sorted(vals, reverse=True)
    if len(vals) <= nth_largest:
        return 0.0
    self_val = int(self_map.get(node, {}).get(metric, 0))
    denom = max(self_val, 1)
    return vals[nth_largest] / denom


def motif_component_passes(subgraph, category):
    self_edges, pair_edges = _component_self_and_pair_edges(subgraph)
    self_map = _self_metric_map(self_edges)

    if category == "self_only":
        return True

    if category == "pair_only":
        if len(pair_edges) != 1:
            return False
        edge = pair_edges[0]
        nodes = [edge["u"], edge["v"]]
        per_node = []
        for node in nodes:
            self_r = int(self_map.get(node, {}).get("support_reads", 0))
            self_u = int(self_map.get(node, {}).get("support_umis", 0))
            read_frac = edge["support_reads"] / (edge["support_reads"] + self_r) if (edge["support_reads"] + self_r) else 0.0
            umi_frac = edge["support_umis"] / (edge["support_umis"] + self_u) if (edge["support_umis"] + self_u) else 0.0
            per_node.append(max(read_frac, umi_frac))
        return min(per_node) >= 0.35 and edge["support_umis"] >= 2

    if category == "triangle_only":
        if len(pair_edges) != 3:
            return False
        if max(
            _pairwise_fraction(pair_edges, self_edges, "support_reads"),
            _pairwise_fraction(pair_edges, self_edges, "support_umis"),
        ) < 0.55:
            return False
        if max(
            _weakest_to_median_ratio(pair_edges, "support_reads"),
            _weakest_to_median_ratio(pair_edges, "support_umis"),
        ) < 0.25:
            return False
        for node in subgraph.nodes():
            if max(
                _node_pair_fraction(subgraph, node, self_map, "support_reads"),
                _node_pair_fraction(subgraph, node, self_map, "support_umis"),
            ) < 0.45:
                return False
        return True

    if category == "clique4_only":
        if len(pair_edges) != 6:
            return False
        if max(
            _pairwise_fraction(pair_edges, self_edges, "support_reads"),
            _pairwise_fraction(pair_edges, self_edges, "support_umis"),
        ) < 0.65:
            return False
        if max(
            _weakest_to_median_ratio(pair_edges, "support_reads"),
            _weakest_to_median_ratio(pair_edges, "support_umis"),
        ) < 0.18:
            return False
        for node in subgraph.nodes():
            if max(
                _nth_incident_pair_to_self(subgraph, node, self_map, "support_reads", 1),
                _nth_incident_pair_to_self(subgraph, node, self_map, "support_umis", 1),
            ) < 0.25:
                return False
        return True

    if category == "clique5_only":
        if len(pair_edges) != 10:
            return False
        if max(
            _pairwise_fraction(pair_edges, self_edges, "support_reads"),
            _pairwise_fraction(pair_edges, self_edges, "support_umis"),
        ) < 0.72:
            return False
        if max(
            _weakest_to_median_ratio(pair_edges, "support_reads"),
            _weakest_to_median_ratio(pair_edges, "support_umis"),
        ) < 0.15:
            return False
        for node in subgraph.nodes():
            if max(
                _nth_incident_pair_to_self(subgraph, node, self_map, "support_reads", 2),
                _nth_incident_pair_to_self(subgraph, node, self_map, "support_umis", 2),
            ) < 0.18:
                return False
        return True

    return False


def prune_edges_structure_aware(
    pc,
    top1_reads_map,
    top1_umis_map,
    top1_alpha,
    top1_alpha_umi,
    require_pass_both_ends=False,
):
    if len(pc) == 0:
        empty = pc.copy()
        stats = {
            "small_components_seen": 0,
            "small_components_motif_kept": 0,
            "edges_kept_by_motif": 0,
            "edges_kept_by_fallback": 0,
        }
        return empty, empty, stats

    graph, _ = build_graph_from_pair_counts(pc, col5="BC5n", col3="BC3n", wcol="support_reads")
    edge_lookup = {}
    for row in pc.itertuples(index=False):
        u = min(row.BC5n, row.BC3n)
        v = max(row.BC5n, row.BC3n)
        edge_lookup[(u, v)] = {
            "BC5n": u,
            "BC3n": v,
            "support_reads": int(row.support_reads),
            "support_umis": int(getattr(row, "support_umis", 0)),
        }

    keep_pairs = {}
    stats = {
        "small_components_seen": 0,
        "small_components_motif_kept": 0,
        "edges_kept_by_motif": 0,
        "edges_kept_by_fallback": 0,
    }

    for comp_nodes in nx.connected_components(graph):
        subgraph = graph.subgraph(comp_nodes).copy()
        category = component_category(subgraph)
        n_nodes = subgraph.number_of_nodes()
        use_motif = category in {"self_only", "pair_only", "triangle_only", "clique4_only", "clique5_only"}
        if use_motif:
            stats["small_components_seen"] += 1
        if use_motif and motif_component_passes(subgraph, category):
            stats["small_components_motif_kept"] += 1
            for u, v, _ in subgraph.edges(data=True):
                key = (min(u, v), max(u, v))
                keep_pairs[key] = True
                stats["edges_kept_by_motif"] += 1
            continue

        for u, v, data in subgraph.edges(data=True):
            key = (min(u, v), max(u, v))
            if u == v:
                keep_pairs[key] = True
                stats["edges_kept_by_fallback"] += 1
                continue
            w_reads = int(data.get("weight", 0))
            w_umis = int(data.get("support_umis", 0))
            pass_u = _endpoint_pass_with_umi(
                w_reads,
                w_umis,
                top1_reads_map.get(u, 0),
                top1_umis_map.get(u, 0),
                top1_alpha,
                top1_alpha_umi,
            )
            pass_v = _endpoint_pass_with_umi(
                w_reads,
                w_umis,
                top1_reads_map.get(v, 0),
                top1_umis_map.get(v, 0),
                top1_alpha,
                top1_alpha_umi,
            )
            if (pass_u and pass_v) if require_pass_both_ends else (pass_u or pass_v):
                keep_pairs[key] = True
                stats["edges_kept_by_fallback"] += 1

    pc_keep = pd.DataFrame(list(keep_pairs.values()))
    if keep_pairs:
        rows = [edge_lookup[k] for k in sorted(keep_pairs.keys())]
        pc_keep = pd.DataFrame(rows)
    else:
        pc_keep = pc.iloc[0:0].copy()

    keep_keys = set(keep_pairs.keys())
    pc_drop = pc[
        ~pc.apply(lambda r: (min(r["BC5n"], r["BC3n"]), max(r["BC5n"], r["BC3n"])) in keep_keys, axis=1)
    ].copy()
    if len(pc_drop):
        pc_drop = pc_drop.assign(drop_reason="structure_aware_or_alpha_umi_prune")

    return pc_keep, pc_drop, stats

def filter_pairs_three_stage(
    df_reads,
    bc5_col="BC5_20bp",
    bc3_col="BC3_20bp_rc",
    umi3_col="putative_umi",
    umi5_col="putative_umi_5p",
    read_id_col="read_id",
    PAIR_MIN=None,
    auto_pair_min_floor=10,
    auto_pair_min_quantile=0.1,
    TOP1_ALPHA=0.3,
    TOP1_ALPHA_UMI=0.3,
    require_pass_both_ends=False,
    drop_umiA_ratio_gt=0.5,
    keep_pair_support_col=True,
):
    df = df_reads.copy()
    df["_b5"] = df[bc5_col].map(norm_bc) if bc5_col in df.columns else ""
    df["_b3"] = df[bc3_col].map(norm_bc) if bc3_col in df.columns else ""
    df["_umi3"] = df[umi3_col].map(norm_umi) if umi3_col in df.columns else ""
    df["_umi5"] = df[umi5_col].map(norm_umi) if umi5_col in df.columns else ""

    n0 = len(df)
    df = df[(df["_b5"] != "") | (df["_b3"] != "")].copy()
    n_drop_both_empty_bc = int(n0 - len(df))

    if umi3_col in df.columns:
        mask_badA = df["_umi3"].apply(lambda s: a_ratio(s) > drop_umiA_ratio_gt)
        n_drop_badA = int(mask_badA.sum())
        df = df[~mask_badA].copy()
    else:
        n_drop_badA = 0

    mask_paired = (df["_b5"] != "") & (df["_b3"] != "")
    df_paired = df[mask_paired].copy()
    df_single_end = df[~mask_paired].copy()
    if len(df_paired) > 0:
        # Treat barcode pairs as unordered groups from the start:
        # (BC5=A, BC3=B) and (BC5=B, BC3=A) belong to the same pair bucket.
        df_paired["pair_u"] = df_paired[["_b5", "_b3"]].min(axis=1)
        df_paired["pair_v"] = df_paired[["_b5", "_b3"]].max(axis=1)
    else:
        df_paired["pair_u"] = ""
        df_paired["pair_v"] = ""

    pair_counts_all = build_pair_counts_with_umis(df_paired)

    PAIR_MIN, pair_min_mode = resolve_pair_min(
        pair_counts_all,
        pair_min=PAIR_MIN,
        auto_pair_min_floor=auto_pair_min_floor,
        auto_pair_min_quantile=auto_pair_min_quantile,
    )

    pair_counts_min_kept = pair_counts_all[pair_counts_all["support_reads"] >= PAIR_MIN].copy()
    pair_counts_min_drop = pair_counts_all[pair_counts_all["support_reads"] < PAIR_MIN].copy()
    pair_counts_min_drop = pair_counts_min_drop.assign(drop_reason=f"pair_min<{PAIR_MIN}")

    df_paired2 = df_paired.merge(
        pair_counts_all.rename(columns={"BC5n": "pair_u", "BC3n": "pair_v"}),
        on=["pair_u", "pair_v"],
        how="left",
    )
    df_paired_min_kept = df_paired2[df_paired2["support_reads"] >= PAIR_MIN].copy()

    pc = pair_counts_min_kept.copy()
    if len(pc) == 0:
        pair_counts_final = pc.copy()
        pair_counts_top1_drop = pc.copy()
        df_paired_final = df_paired_min_kept.iloc[0:0].copy()
        motif_stats = {
            "small_components_seen": 0,
            "small_components_motif_kept": 0,
            "edges_kept_by_motif": 0,
            "edges_kept_by_fallback": 0,
        }
    else:
        top1_reads_map = _top1_map_from_metric(pc, "support_reads")
        top1_umis_map = _top1_map_from_metric(pc, "support_umis")
        pair_counts_final, pair_counts_top1_drop, motif_stats = prune_edges_structure_aware(
            pc,
            top1_reads_map=top1_reads_map,
            top1_umis_map=top1_umis_map,
            top1_alpha=TOP1_ALPHA,
            top1_alpha_umi=TOP1_ALPHA_UMI,
            require_pass_both_ends=require_pass_both_ends,
        )
        df_paired_final = df_paired_min_kept.merge(
            pair_counts_final.rename(columns={"BC5n": "pair_u", "BC3n": "pair_v"})[["pair_u", "pair_v"]],
            on=["pair_u", "pair_v"],
            how="inner",
        )

    if not keep_pair_support_col and "support_reads" in df_paired_final.columns:
        df_paired_final = df_paired_final.drop(columns=["support_reads"], errors="ignore")

    df_final = pd.concat([df_single_end, df_paired_final], ignore_index=True)
    dropped_pairs = pd.concat([pair_counts_min_drop, pair_counts_top1_drop], ignore_index=True)

    stats = {
        "rows_in": int(n0),
        "drop_both_empty_bc": int(n_drop_both_empty_bc),
        "drop_putative_umi_A_ratio_gt": float(drop_umiA_ratio_gt),
        "drop_putative_umi_badA_rows": int(n_drop_badA),
        "single_end_reads": int(len(df_single_end)),
        "paired_reads_before": int(len(df_paired)),
        "pairs_total": int(len(pair_counts_all)),
        "pairs_kept_pair_min": int(len(pair_counts_min_kept)),
        "pairs_dropped_pair_min": int(len(pair_counts_min_drop)),
        "pairs_kept_final": int(len(pair_counts_final)),
        "pairs_dropped_top1": int(len(pair_counts_top1_drop)),
        "paired_reads_kept_final": int(len(df_paired_final)),
        "rows_out_final": int(len(df_final)),
        "require_pass_both_ends": bool(require_pass_both_ends),
        "PAIR_MIN": int(PAIR_MIN),
        "PAIR_MIN_mode": pair_min_mode,
        "auto_pair_min_floor": int(auto_pair_min_floor),
        "auto_pair_min_quantile": float(auto_pair_min_quantile),
        "TOP1_ALPHA": float(TOP1_ALPHA),
        "TOP1_ALPHA_UMI": float(TOP1_ALPHA_UMI),
        "structure_small_components_seen": int(motif_stats["small_components_seen"]),
        "structure_small_components_motif_kept": int(motif_stats["small_components_motif_kept"]),
        "structure_edges_kept_by_motif": int(motif_stats["edges_kept_by_motif"]),
        "structure_edges_kept_by_fallback": int(motif_stats["edges_kept_by_fallback"]),
    }

    for dfx in (df_final, df_single_end, df_paired_final):
        dfx.drop(columns=["_b5", "_b3", "_umi3", "_umi5", "pair_u", "pair_v"], errors="ignore", inplace=True)

    return (
        df_final,
        df_single_end,
        df_paired_final,
        pair_counts_all,
        pair_counts_min_kept,
        pair_counts_final,
        dropped_pairs,
        stats,
    )
