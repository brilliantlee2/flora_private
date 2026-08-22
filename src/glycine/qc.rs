use super::utils::ProcessedRecord;
use fxread::Record;
use std::sync::OnceLock;

static QUALITY_ERROR_PROBABILITIES: OnceLock<[f64; 256]> = OnceLock::new();

fn quality_error_probabilities() -> &'static [f64; 256] {
    QUALITY_ERROR_PROBABILITIES.get_or_init(|| {
        std::array::from_fn(|quality| 10.0_f64.powf((quality as f64 - 33.0_f64) / (-10.0_f64)))
    })
}

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

    let error_probabilities = quality_error_probabilities();
    let read_err_rate: f64 = raw_record
        .qual()
        .unwrap()
        .iter()
        .map(|&v| error_probabilities[v as usize])
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

#[cfg(test)]
mod tests {
    use super::quality_error_probabilities;

    #[test]
    fn quality_lookup_matches_original_powf_formula() {
        let lookup = quality_error_probabilities();
        for quality in 0_u16..=255 {
            let expected = 10.0_f64.powf((quality as f64 - 33.0_f64) / (-10.0_f64));
            assert_eq!(lookup[quality as usize], expected);
        }
    }
}
