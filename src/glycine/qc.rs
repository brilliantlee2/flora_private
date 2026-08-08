use super::utils::ProcessedRecord;
use fxread::Record;

// 根据长度和Q值对reads做QC
pub fn qc(
    raw_record: &Record,
    raw_record_len: usize,
    min_len: usize,
    min_qual: f64,
) -> Option<(ProcessedRecord, usize)> {
    if raw_record_len < min_len {
        let length_filtered_record_id: &Vec<u8> =
            &[b"@", raw_record.id(), b" short-length"].concat();
        let length_filtered_record = ProcessedRecord::new(
            length_filtered_record_id.to_owned(),
            raw_record.seq().to_vec(),
            raw_record.plus().unwrap().to_vec(),
            raw_record.qual().unwrap().to_vec(),
        );

        return Some((length_filtered_record, 16));
    }

    let read_err_rate: f64 = raw_record
        .qual()
        .unwrap()
        .iter()
        .map(|&v| 10.0_f64.powf((v as f64 - 33.0_f64) / (-10.0_f64)))
        .sum::<f64>()
        / raw_record_len as f64;
    let read_qval: f64 = -10.0_f64 * read_err_rate.log10();

    if read_qval < min_qual {
        let qval_filtered_record_id: &Vec<u8> = &[b"@", raw_record.id(), b" low-quality"].concat();
        let qval_filtered_record = ProcessedRecord::new(
            qval_filtered_record_id.to_owned(),
            raw_record.seq().to_vec(),
            raw_record.plus().unwrap().to_vec(),
            raw_record.qual().unwrap().to_vec(),
        );

        return Some((qval_filtered_record, 17));
    }

    None
}
