mod args;
mod file_system;
mod identifier;
mod qc;
mod reader;
mod utils;
mod writer;

use args::*;
use file_system::*;
use flate2::write::GzEncoder;
use fxread::Record;
use identifier::*;
use reader::*;
use std::fs::File;
use std::path::PathBuf;
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;
use utils::*;
use writer::*;

pub fn run_from<I, T>(args: I) -> anyhow::Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let matches = parse_argument_from(args);
    let fq_file_path = matches.get_one::<String>("fastq file").unwrap();
    let tso_seq: Vec<u8> = matches
        .get_one::<String>("tso seq")
        .unwrap()
        .bytes()
        .collect();
    let rtp_seq: Vec<u8> = matches
        .get_one::<String>("rtp seq")
        .unwrap()
        .bytes()
        .collect();
    let output_dir_path = matches.get_one::<PathBuf>("output dir").unwrap();
    let sample = matches.get_one::<String>("sample").unwrap();
    let err_threshold_vec: Vec<ErrThreshold> = matches
        .get_one::<String>("err threshold")
        .unwrap()
        .split(',')
        .map(|s| parse_number(s))
        .collect();
    let shift_threshold_vec: Vec<usize> = matches
        .get_one::<String>("shift threshold")
        .unwrap()
        .split(',')
        .map(|s| s.parse::<usize>().unwrap())
        .collect();
    let min_len = *matches.get_one::<usize>("min len").unwrap();
    let min_qual = *matches.get_one::<f64>("min qual").unwrap();
    let trim_len = *matches.get_one::<usize>("trim len").unwrap();
    let tail_len = *matches.get_one::<usize>("tail len").unwrap();
    let umi_len = *matches.get_one::<usize>("umi len").unwrap();
    let thread_num = *matches.get_one::<usize>("thread").unwrap();
    let keep_all_outputs = matches.get_flag("keep all outputs");

    if trim_len >= tso_seq.len().min(rtp_seq.len()) {
        panic!("Trim_len must be less than TSO and RTP seq. Please adjust the parameters!");
    }

    let min_len = tso_seq.len().max(rtp_seq.len()).max(tail_len).max(min_len);

    let fq_reader = read_fq_file(fq_file_path).unwrap();
    let mut encoders: Vec<GzEncoder<File>> =
        create_output_dir(output_dir_path, sample, keep_all_outputs).unwrap();

    let pool = ThreadPool::new(thread_num);

    let tso_seq: Arc<Vec<u8>> = Arc::new(tso_seq);
    let rtp_seq: Arc<Vec<u8>> = Arc::new(rtp_seq);

    let (err_threshold_start, err_threshold_end): (ErrThreshold, ErrThreshold);
    if err_threshold_vec.len() < 2 {
        err_threshold_start = err_threshold_vec[0].clone();
        err_threshold_end = err_threshold_vec[0].clone();
    } else {
        err_threshold_start = err_threshold_vec[0].clone();
        err_threshold_end = err_threshold_vec[1].clone();
    }

    let err_threshold_tso = convert_ld(&err_threshold_start, tso_seq.len());
    let err_threshold_rtp = convert_ld(&err_threshold_start, rtp_seq.len());
    let err_threshold_tso_comp = convert_ld(&err_threshold_end, tso_seq.len() - trim_len);
    let err_threshold_rtp_comp = convert_ld(&err_threshold_end, rtp_seq.len() - trim_len);

    let (shift_threshold_start, shift_threshold_end): (usize, usize);
    if shift_threshold_vec.len() < 2 {
        shift_threshold_start = shift_threshold_vec[0];
        shift_threshold_end = shift_threshold_vec[0];
    } else {
        shift_threshold_start = shift_threshold_vec[0];
        shift_threshold_end = shift_threshold_vec[1];
    }

    let mut summary_read_num: Vec<u64> = vec![0; 5];
    let mut non_chimeric_read_num: Vec<u64> = vec![0; 16];
    let mut chimeric_read_num: Vec<u64> = vec![0; 3];
    let mut total_read_num: u64 = 0;
    let mut total_base_num: u64 = 0;
    let mut valid_base_num: u64 = 0;
    let mut rescued_read_num: u64 = 0;
    let mut full_length_plus_rescued_num: u64 = 0;

    let mut batch: Vec<Record> = vec![];

    for raw_record in fq_reader {
        total_read_num += 1;
        total_base_num += raw_record.seq().len() as u64;

        batch.push(raw_record);

        if total_read_num % 100_000 == 0 {
            let (tx, rx) = mpsc::channel();

            for raw_record in batch {
                let tx = tx.clone();
                let tso_seq = Arc::clone(&tso_seq);
                let rtp_seq = Arc::clone(&rtp_seq);

                pool.execute(move || {
                    let res_fq_result = classify_reads(
                        raw_record,
                        &tso_seq,
                        &rtp_seq,
                        err_threshold_tso,
                        err_threshold_rtp,
                        err_threshold_tso_comp,
                        err_threshold_rtp_comp,
                        shift_threshold_start,
                        shift_threshold_end,
                        min_len,
                        min_qual,
                        trim_len,
                        tail_len,
                        umi_len,
                    );
                    tx.send(res_fq_result).unwrap_or_else(|e| {
                        eprintln!("Failed to send data on the channel: {}", e);
                    });
                });
            }

            drop(tx);

            for (res_fq_vec, base_num) in rx.into_iter() {
                let first_idx = res_fq_vec[0].1;

                if first_idx <= 15 {
                    summary_read_num[3] += 1;
                    if first_idx <= 3 {
                        summary_read_num[4] += 1;
                    }
                } else if first_idx == 16 {
                    summary_read_num[0] += 1;
                } else if first_idx == 17 {
                    summary_read_num[1] += 1;
                } else if first_idx == 18 {
                    summary_read_num[2] += 1;
                } else {
                    summary_read_num[2] += 1;
                    summary_read_num[4] += 1;
                    rescued_read_num += 1;
                }

                for (res_fq_content, res_fq_idx) in res_fq_vec.into_iter() {
                    let file_path_idx: usize;

                    if res_fq_idx <= 3 {
                        non_chimeric_read_num[res_fq_idx] += 1;
                        valid_base_num += base_num as u64;
                        file_path_idx = 0;
                    } else if res_fq_idx <= 12 {
                        non_chimeric_read_num[res_fq_idx] += 1;
                        file_path_idx = 1;
                    } else if res_fq_idx <= 15 {
                        non_chimeric_read_num[res_fq_idx] += 1;
                        file_path_idx = 2;
                    } else if res_fq_idx <= 17 {
                        file_path_idx = 3;
                    } else if res_fq_idx == 18 {
                        chimeric_read_num[2] += 1;
                        file_path_idx = 2;
                    } else if res_fq_idx == 19 {
                        chimeric_read_num[0] += 1;
                        valid_base_num += base_num as u64;
                        file_path_idx = 4;
                    } else {
                        chimeric_read_num[1] += 1;
                        valid_base_num += base_num as u64;
                        file_path_idx = 4;
                    }

                    let merge_into_full_length_plus_rescued =
                        res_fq_idx <= 3 || res_fq_idx == 19 || res_fq_idx == 20;
                    let merged_record = if merge_into_full_length_plus_rescued {
                        Some(res_fq_content.clone())
                    } else {
                        None
                    };

                    if keep_all_outputs {
                        write_fq_file(&mut encoders[file_path_idx], res_fq_content).unwrap();
                    }
                    if let Some(merged_record) = merged_record {
                        full_length_plus_rescued_num += 1;
                        let merged_idx = if keep_all_outputs { 5 } else { 0 };
                        write_fq_file(&mut encoders[merged_idx], merged_record).unwrap();
                    }
                }
            }

            pool.join();

            flush_fq_file(&mut encoders).unwrap();

            batch = vec![];
        }
    }

    if !batch.is_empty() {
        let (tx, rx) = mpsc::channel();

        for raw_record in batch {
            let tx = tx.clone();
            let tso_seq = Arc::clone(&tso_seq);
            let rtp_seq = Arc::clone(&rtp_seq);

            pool.execute(move || {
                let res_fq_result = classify_reads(
                    raw_record,
                    &tso_seq,
                    &rtp_seq,
                    err_threshold_tso,
                    err_threshold_rtp,
                    err_threshold_tso_comp,
                    err_threshold_rtp_comp,
                    shift_threshold_start,
                    shift_threshold_end,
                    min_len,
                    min_qual,
                    trim_len,
                    tail_len,
                    umi_len,
                );
                tx.send(res_fq_result).unwrap_or_else(|e| {
                    eprintln!("Failed to send data on the channel: {}", e);
                });
            });
        }

        drop(tx);

        for (res_fq_vec, base_num) in rx.into_iter() {
            let first_idx = res_fq_vec[0].1;

            if first_idx <= 15 {
                summary_read_num[3] += 1;
                if first_idx <= 3 {
                    summary_read_num[4] += 1;
                }
            } else if first_idx == 16 {
                summary_read_num[0] += 1;
            } else if first_idx == 17 {
                summary_read_num[1] += 1;
            } else if first_idx == 18 {
                summary_read_num[2] += 1;
            } else {
                summary_read_num[2] += 1;
                summary_read_num[4] += 1;
                rescued_read_num += 1;
            }

            for (res_fq_content, res_fq_idx) in res_fq_vec.into_iter() {
                let file_path_idx: usize;

                if res_fq_idx <= 3 {
                    non_chimeric_read_num[res_fq_idx] += 1;
                    valid_base_num += base_num as u64;
                    file_path_idx = 0;
                } else if res_fq_idx <= 12 {
                    non_chimeric_read_num[res_fq_idx] += 1;
                    file_path_idx = 1;
                } else if res_fq_idx <= 15 {
                    non_chimeric_read_num[res_fq_idx] += 1;
                    file_path_idx = 2;
                } else if res_fq_idx <= 17 {
                    file_path_idx = 3;
                } else if res_fq_idx == 18 {
                    chimeric_read_num[2] += 1;
                    file_path_idx = 2;
                } else if res_fq_idx == 19 {
                    chimeric_read_num[0] += 1;
                    valid_base_num += base_num as u64;
                    file_path_idx = 4;
                } else {
                    chimeric_read_num[1] += 1;
                    valid_base_num += base_num as u64;
                    file_path_idx = 4;
                }

                let merge_into_full_length_plus_rescued =
                    res_fq_idx <= 3 || res_fq_idx == 19 || res_fq_idx == 20;
                let merged_record = if merge_into_full_length_plus_rescued {
                    Some(res_fq_content.clone())
                } else {
                    None
                };

                if keep_all_outputs {
                    write_fq_file(&mut encoders[file_path_idx], res_fq_content).unwrap();
                }
                if let Some(merged_record) = merged_record {
                    full_length_plus_rescued_num += 1;
                    let merged_idx = if keep_all_outputs { 5 } else { 0 };
                    write_fq_file(&mut encoders[merged_idx], merged_record).unwrap();
                }
            }
        }

        pool.join();

        flush_fq_file(&mut encoders).unwrap();
    }

    summary_read_num.insert(0, total_read_num);
    let full_length_num: u64 = non_chimeric_read_num[..4].iter().sum();
    let non_full_length_num: u64 = non_chimeric_read_num[4..13].iter().sum();
    let discarded_num: u64 = non_chimeric_read_num[13..].iter().sum();
    non_chimeric_read_num.insert(0, full_length_num);
    non_chimeric_read_num.insert(5, non_full_length_num);
    non_chimeric_read_num.insert(15, discarded_num);

    let mut statistic_content = format!("Summary\nTotal_base_count\tValid_base_count\tValid_base_proportion(%)\n{}\t{}\t{:.2}\nType\tRead_count\tRead_proportion(%)", total_base_num, valid_base_num, valid_base_num as f64 * 100.0 / total_base_num as f64);

    let summary_type_vec: Vec<&str> = vec![
        "Total",
        "Length-filtered",
        "QC-filtered",
        "Full-length+rescued",
        "Chimeric",
        "Non-chimeric",
        "Full-length",
    ];
    let summary_count_vec: Vec<u64> = vec![
        total_read_num,
        summary_read_num[1],
        summary_read_num[2],
        full_length_plus_rescued_num,
        summary_read_num[3],
        summary_read_num[4],
        summary_read_num[5],
    ];

    for (t, &n) in summary_type_vec.into_iter().zip(summary_count_vec.iter()) {
        statistic_content = format!(
            "{}\n{}\t{}\t{:.2}",
            statistic_content,
            t,
            n,
            n as f64 * 100.0 / total_read_num as f64
        );
    }

    statistic_content = format!(
        "{}\n\nNon-chimeric\nType\tRead_count\tRead_proportion(%)\nTotal\t{}\t100.00",
        statistic_content, summary_read_num[4]
    );

    let non_chimeric_type_vec: Vec<&str> = vec![
        "Full-length",
        "High-confidence strand:+",
        "High-confidence strand:-",
        "Low-confidence strand:+",
        "Low-confidence strand:-",
        "Non-full-length",
        "Double-rtp-double-polya/t",
        "Non-tso-polya-rtp",
        "Rtp-polyt-non-tso",
        "Tso-rtp-non-polya/t",
        "Rtp-tso-non-polya/t",
        "Single-5'-tso",
        "Single-5'-rtp",
        "Single-3'-tso",
        "Single-3'-rtp",
        "Discarded",
        "Double-tso",
        "Double-rtp",
        "Non-primer",
    ];

    for (t, n) in non_chimeric_type_vec
        .into_iter()
        .zip(non_chimeric_read_num.into_iter())
    {
        statistic_content = format!(
            "{}\n{}\t{}\t{:.2}",
            statistic_content,
            t,
            n,
            n as f64 * 100.0 / summary_read_num[4] as f64
        );
    }

    statistic_content = format!(
        "{}\n\nChimeric\nChimeric_read_count\tRescued_read_count\tRescued_read_proportion(%)\n{}\t{}\t{:.2}\nType\tRead_count\tRead_proportion(%)\nTotal\t{}\t100.00",
        statistic_content, summary_read_num[3], rescued_read_num, rescued_read_num as f64 * 100.0 / summary_read_num[3] as f64, chimeric_read_num.iter().sum::<u64>()
    );

    let chimeric_type_vec: Vec<&str> = vec!["Rescued strand:+", "Rescued strand:-", "Irrescuable"];

    for (t, &n) in chimeric_type_vec.into_iter().zip(chimeric_read_num.iter()) {
        statistic_content = format!(
            "{}\n{}\t{}\t{:.2}",
            statistic_content,
            t,
            n,
            n as f64 * 100.0 / chimeric_read_num.iter().sum::<u64>() as f64
        );
    }

    statistic_content.push_str("\n");
    create_and_write_statistic_file(output_dir_path, sample, statistic_content)?;
    Ok(())
}
