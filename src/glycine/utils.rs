use edlib_rs::edlibrs::{
    edlibAlignRs, EdlibAlignConfigRs, EdlibAlignModeRs, EdlibAlignTaskRs, EDLIB_STATUS_OK,
};

// 将序列转换为反向序列
pub fn convert_rev(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| c).collect()
}

// 将序列转换为互补序列
pub fn convert_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|c| if c & 2 == 0 { c ^ 21 } else { c ^ 4 })
        .collect()
}

// 将序列转换为反向互补序列
pub fn convert_rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| if c & 2 == 0 { c ^ 21 } else { c ^ 4 })
        .collect()
}

// 生成切除tso/rtp/polyA/polyT后的record
// #[inline]
// fn generate_record(id: &[u8], seq: &[u8], plus: &[u8], qual: &[u8]) -> Vec<u8> {
//     vec![
//         &[&[b'@'], id].concat(),
//         seq,
//         plus,
//         &[qual, &[b'\n']].concat(),
//     ]
//     .join(&b'\n')
// }

// 允许ErrThreshold为usize或f64，分别代表LD阈值和测序错误率阈值
#[derive(Clone)]
pub enum ErrThreshold {
    LD(usize),
    ErrRate(f64),
}

// 解析ErrThreshold为usize或f64
pub fn parse_number(s: &str) -> ErrThreshold {
    if let Ok(i) = s.parse::<usize>() {
        ErrThreshold::LD(i)
    } else if let Ok(f) = s.parse::<f64>() {
        ErrThreshold::ErrRate(f)
    } else {
        panic!("Not a valid usize or f64");
    }
}

// 将ErrThreshold转换为LD
pub fn convert_ld(err_threshold: &ErrThreshold, seq_len: usize) -> usize {
    let err_threshold_res = match err_threshold {
        ErrThreshold::LD(i) => *i,
        ErrThreshold::ErrRate(f) => (seq_len as f64 * f).round() as usize,
    };

    err_threshold_res
}

// 计算最小编辑距离和坐标
pub fn min_edit_distance_and_coord(pattern: &[u8], text: &[u8]) -> (usize, usize) {
    if pattern.is_empty() {
        return (0, text.len());
    }
    if text.is_empty() {
        return (pattern.len(), 0);
    }

    let config = EdlibAlignConfigRs::new(
        -1,
        EdlibAlignModeRs::EDLIB_MODE_HW,
        EdlibAlignTaskRs::EDLIB_TASK_LOC,
        &[],
    );
    let result = edlibAlignRs(pattern, text, &config);
    assert_eq!(result.status, EDLIB_STATUS_OK, "Edlib HW alignment failed");
    assert!(result.editDistance >= 0, "Edlib did not return a distance");

    let end_exclusive = result
        .endLocations
        .as_ref()
        .and_then(|locations| locations.iter().copied().max())
        .map(|location| {
            if location < 0 {
                0
            } else {
                location as usize + 1
            }
        })
        .expect("Edlib did not return an end location");

    (result.editDistance as usize, end_exclusive)
}

// 定义处理后的record结构体
#[derive(Clone)]
pub struct ProcessedRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub plus: Vec<u8>,
    pub qual: Vec<u8>,
}

impl ProcessedRecord {
    // 实现从record信息创建ProcessedRecord的方法
    pub fn new(id: Vec<u8>, seq: Vec<u8>, plus: Vec<u8>, qual: Vec<u8>) -> Self {
        ProcessedRecord {
            id,
            seq,
            plus,
            qual,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::min_edit_distance_and_coord;
    use ndarray::prelude::*;

    fn reference_min_edit_distance_and_coord(pattern: &[u8], text: &[u8]) -> (usize, usize) {
        let pattern_len = pattern.len();
        let text_len = text.len();
        let mut matrix = Array::<usize, _>::zeros((pattern_len + 1, text_len + 1));

        for i in 1..=pattern_len {
            matrix[[i, 0]] = i;
        }
        for i in 1..=pattern_len {
            for j in 1..=text_len {
                if pattern[i - 1] == text[j - 1] {
                    matrix[[i, j]] = matrix[[i - 1, j - 1]];
                } else {
                    matrix[[i, j]] = matrix[[i - 1, j - 1]]
                        .min(matrix[[i - 1, j]])
                        .min(matrix[[i, j - 1]])
                        + 1;
                }
            }
        }

        let last_row = matrix.index_axis(Axis(0), pattern_len);
        let mut min_ed = *last_row.first().unwrap();
        let mut min_idx = 0;
        for (idx, &ed) in last_row.into_iter().enumerate() {
            if ed <= min_ed {
                min_ed = ed;
                min_idx = idx;
            }
        }
        (min_ed, min_idx)
    }

    fn assert_matches_reference(pattern: &[u8], text: &[u8]) {
        assert_eq!(
            min_edit_distance_and_coord(pattern, text),
            reference_min_edit_distance_and_coord(pattern, text),
            "pattern={} text={}",
            String::from_utf8_lossy(pattern),
            String::from_utf8_lossy(text)
        );
    }

    #[test]
    fn edlib_matches_reference_for_required_edge_cases() {
        assert_matches_reference(b"ACGT", b"TTACGTGG");
        assert_matches_reference(b"ACGA", b"TTACGTGG");
        assert_matches_reference(b"ACGTT", b"TTACGTGG");
        assert_matches_reference(b"ACT", b"TTACGTGG");
        assert_matches_reference(b"A", b"AAA");
        assert_eq!(min_edit_distance_and_coord(b"A", b"AAA"), (0, 3));
        assert_matches_reference(b"", b"ACGT");
        assert_matches_reference(b"ACGT", b"");
    }

    #[test]
    fn edlib_matches_reference_for_random_dna() {
        let alphabet = [b'A', b'C', b'G', b'T'];
        let mut state = 0x4d59_5df4_d0f3_3173_u64;
        let mut next_base = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            alphabet[((state >> 32) as usize) & 3]
        };

        for case_idx in 0..2_000 {
            let pattern_len = case_idx % 35;
            let text_len = (case_idx * 17) % 81;
            let pattern: Vec<u8> = (0..pattern_len).map(|_| next_base()).collect();
            let text: Vec<u8> = (0..text_len).map(|_| next_base()).collect();
            assert_matches_reference(&pattern, &text);
        }
    }
}
