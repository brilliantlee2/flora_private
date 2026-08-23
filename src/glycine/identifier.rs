use super::qc::*;
use super::utils::*;
use edlib_rs::edlibrs::*;
use fxread::Record;
use std::ops::Range;

fn chimeric_tail_search_span(tail_len: usize, umi_len: usize) -> usize {
    (tail_len * 3).max(umi_len + tail_len + 10)
}

fn build_rescued_record_id(raw_id: &[u8], record_number: &[u8], strand: &[u8]) -> Vec<u8> {
    [
        b"@",
        raw_id,
        b"_rescued_No",
        record_number,
        b"_strand_",
        strand,
    ]
    .concat()
}

fn checked_rescue_range(start: usize, end: usize, read_len: usize) -> Option<Range<usize>> {
    (start < end && end <= read_len).then_some(start..end)
}

fn plus_rescue_range(
    tso_end: usize,
    rtp_rev_comp_start: usize,
    polya_offset: usize,
    read_len: usize,
) -> Option<Range<usize>> {
    let start = tso_end.checked_add(1)?;
    let end = rtp_rev_comp_start.checked_sub(polya_offset)?;
    checked_rescue_range(start, end, read_len)
}

fn minus_rescue_range(
    rtp_end: usize,
    polyt_offset: usize,
    tso_rev_comp_start: usize,
    read_len: usize,
) -> Option<Range<usize>> {
    let start = rtp_end.checked_add(polyt_offset)?.checked_add(1)?;
    checked_rescue_range(start, tso_rev_comp_start, read_len)
}

pub struct PrimerSequences {
    tso: Vec<u8>,
    rtp: Vec<u8>,
    tso_rev_comp: Vec<u8>,
    rtp_rev_comp: Vec<u8>,
    tso_trimmed_comp: Vec<u8>,
    rtp_trimmed_comp: Vec<u8>,
}

impl PrimerSequences {
    pub fn new(tso: Vec<u8>, rtp: Vec<u8>, trim_len: usize) -> Self {
        let tso_rev_comp = convert_rev_comp(&tso);
        let rtp_rev_comp = convert_rev_comp(&rtp);
        let tso_trimmed_comp = convert_comp(&tso[trim_len..]);
        let rtp_trimmed_comp = convert_comp(&rtp[trim_len..]);
        Self {
            tso,
            rtp,
            tso_rev_comp,
            rtp_rev_comp,
            tso_trimmed_comp,
            rtp_trimmed_comp,
        }
    }
}

// 识别全长转录本、非全长转录本和discarded等
pub fn classify_reads(
    raw_record: Record,
    primers: &PrimerSequences,
    err_threshold_tso: usize,
    err_threshold_rtp: usize,
    err_threshold_tso_comp: usize,
    err_threshold_rtp_comp: usize,
    shift_threshold_start: usize,
    shift_threshold_end: usize,
    min_len: usize,
    min_qual: f64,
    tail_len: usize,
    umi_len: usize,
) -> (Vec<(ProcessedRecord, usize)>, usize) {
    let tso_seq = &primers.tso;
    let rtp_seq = &primers.rtp;
    let raw_record_len = raw_record.seq().len();
    let tail_search_span = chimeric_tail_search_span(tail_len, umi_len);

    if let Some((discarded_record, discarded_idx)) =
        qc(&raw_record, raw_record_len, min_len, min_qual)
    {
        return (vec![(discarded_record, discarded_idx)], 0);
    }

    let (mut tso_config, mut rtp_config) =
        (EdlibAlignConfigRs::default(), EdlibAlignConfigRs::default());
    tso_config.k = err_threshold_tso as i32;
    tso_config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    tso_config.task = EdlibAlignTaskRs::EDLIB_TASK_LOC;
    rtp_config.k = err_threshold_rtp as i32;
    rtp_config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    rtp_config.task = EdlibAlignTaskRs::EDLIB_TASK_LOC;

    let tso_align_res = edlibAlignRs(tso_seq, raw_record.seq(), &tso_config);
    let tso_rev_comp_align_res = edlibAlignRs(&primers.tso_rev_comp, raw_record.seq(), &tso_config);
    let rtp_align_res = edlibAlignRs(rtp_seq, raw_record.seq(), &rtp_config);
    let rtp_rev_comp_align_res = edlibAlignRs(&primers.rtp_rev_comp, raw_record.seq(), &rtp_config);

    let mut primer_loc_vec: Vec<(usize, usize)> = vec![];

    if tso_align_res.endLocations.is_some() {
        for loc in tso_align_res.endLocations.unwrap() {
            primer_loc_vec.push((if loc < 0 { 0 } else { loc as usize }, 0));
        }
    }

    if rtp_align_res.endLocations.is_some() {
        for loc in rtp_align_res.endLocations.unwrap() {
            primer_loc_vec.push((if loc < 0 { 0 } else { loc as usize }, 1));
        }
    }

    if rtp_rev_comp_align_res.startLocations.is_some() {
        for loc in rtp_rev_comp_align_res.startLocations.unwrap() {
            primer_loc_vec.push((if loc < 0 { 0 } else { loc as usize }, 2));
        }
    }

    if tso_rev_comp_align_res.startLocations.is_some() {
        for loc in tso_rev_comp_align_res.startLocations.unwrap() {
            primer_loc_vec.push((if loc < 0 { 0 } else { loc as usize }, 3));
        }
    }

    if primer_loc_vec.is_empty() {
        let non_primer_record_id: Vec<u8> = [b"@", raw_record.id(), b" non-primer"].concat();
        let non_primer_record = ProcessedRecord::new(
            non_primer_record_id,
            raw_record.seq().to_vec(),
            raw_record.plus().unwrap().to_vec(),
            raw_record.qual().unwrap().to_vec(),
        );

        return (vec![(non_primer_record, 15)], 0);
    }

    let mut primer_loc_cluster_vec = vec![];
    let primer_len_diff_threshold =
        tso_seq.len().max(rtp_seq.len()) - tso_seq.len().min(rtp_seq.len()) + 5;

    primer_loc_vec.sort_by(|a, b| a.0.cmp(&b.0));

    for (loc, id) in primer_loc_vec.into_iter() {
        if primer_loc_cluster_vec.is_empty() {
            primer_loc_cluster_vec.push(vec![(loc, id)]);
        }

        let last_loc = primer_loc_cluster_vec.last().unwrap().last().unwrap().0;
        if loc - last_loc <= primer_len_diff_threshold {
            primer_loc_cluster_vec.last_mut().unwrap().push((loc, id));
        } else {
            primer_loc_cluster_vec.push(vec![(loc, id)]);
        }
    }

    if primer_loc_cluster_vec.len() > 2
        || (primer_loc_cluster_vec.len() == 2
            && primer_loc_cluster_vec[0][0].0 >= shift_threshold_start
            && primer_loc_cluster_vec[1].last().unwrap().0 + shift_threshold_end < raw_record_len)
    {
        let mut fq_content: Vec<(ProcessedRecord, usize)> = vec![];
        let mut rescued_base_num: usize = 0;

        for start_end_vec in primer_loc_cluster_vec.windows(2) {
            let (
                mut tso_loc_option,
                mut rtp_loc_option,
                mut rtp_rev_comp_loc_option,
                mut tso_rev_comp_loc_option,
            ) = (None, None, None, None);

            for &(start_loc, start_id) in &start_end_vec[0] {
                if start_id == 0 {
                    tso_loc_option = Some(start_loc);
                } else if start_id == 1 {
                    rtp_loc_option = Some(start_loc);
                }
            }

            for &(end_loc, end_id) in &start_end_vec[1] {
                if end_id == 2 {
                    rtp_rev_comp_loc_option = Some(end_loc);
                } else if end_id == 3 {
                    tso_rev_comp_loc_option = Some(end_loc);
                }
            }

            if rtp_loc_option.is_some() && rtp_rev_comp_loc_option.is_some() {
                if rtp_loc_option.unwrap() + 1 < raw_record_len
                    && rtp_rev_comp_loc_option.unwrap() > 0
                {
                    let rtp_downstream_seq =
                        if rtp_loc_option.unwrap() + tail_search_span + 1 < raw_record_len {
                            raw_record.seq()[rtp_loc_option.unwrap() + 1
                                ..=rtp_loc_option.unwrap() + tail_search_span]
                                .to_vec()
                        } else {
                            raw_record.seq()[rtp_loc_option.unwrap() + 1..].to_vec()
                        };

                    let polyt_loc_option = rtp_downstream_seq
                        .windows(tail_len)
                        .enumerate()
                        .find(|(_, window)| window.iter().all(|&b| b == 84))
                        .map(|(idx, _)| idx);

                    let rtp_rev_comp_upstream_seq =
                        if rtp_rev_comp_loc_option.unwrap() > tail_search_span {
                            raw_record.seq()[rtp_rev_comp_loc_option.unwrap() - tail_search_span
                                ..rtp_rev_comp_loc_option.unwrap()]
                                .to_vec()
                        } else {
                            raw_record.seq()[..rtp_rev_comp_loc_option.unwrap()].to_vec()
                        };

                    let polya_loc_option = rtp_rev_comp_upstream_seq
                        .windows(tail_len)
                        .rev()
                        .enumerate()
                        .find(|(_, window)| window.iter().all(|&b| b == 65))
                        .map(|(idx, _)| idx);

                    if polyt_loc_option.is_some() && polya_loc_option.is_some() {
                        continue;
                    }
                }
            }

            if tso_loc_option.is_some() && rtp_rev_comp_loc_option.is_some() {
                if rtp_rev_comp_loc_option.unwrap() > 0 {
                    let rtp_rev_comp_upstream_seq =
                        if rtp_rev_comp_loc_option.unwrap() > tail_search_span {
                            raw_record.seq()[rtp_rev_comp_loc_option.unwrap() - tail_search_span
                                ..rtp_rev_comp_loc_option.unwrap()]
                                .to_vec()
                        } else {
                            raw_record.seq()[..rtp_rev_comp_loc_option.unwrap()].to_vec()
                        };

                    let polya_loc_option = rtp_rev_comp_upstream_seq
                        .windows(tail_len)
                        .rev()
                        .enumerate()
                        .find(|(_, window)| window.iter().all(|&b| b == 65))
                        .map(|(idx, _)| idx);
                    if polya_loc_option.is_some() {
                        let polya_loc = polya_loc_option.unwrap();
                        if let Some(range) = plus_rescue_range(
                            tso_loc_option.unwrap(),
                            rtp_rev_comp_loc_option.unwrap(),
                            polya_loc,
                            raw_record_len,
                        ) {
                            let record_number = (fq_content.len() + 1)
                                .to_string()
                                .chars()
                                .map(|c| c as u8)
                                .collect::<Vec<u8>>();

                            let rescued_record_id =
                                build_rescued_record_id(raw_record.id(), &record_number, b"plus");
                            let rescued_record = ProcessedRecord::new(
                                rescued_record_id,
                                raw_record.seq()[range.clone()].to_vec(),
                                raw_record.plus().unwrap().to_vec(),
                                raw_record.qual().unwrap()[range].to_vec(),
                            );

                            fq_content.push((rescued_record, 19));
                            rescued_base_num += rtp_rev_comp_loc_option.unwrap()
                                - tso_loc_option.unwrap()
                                + tso_seq.len()
                                + rtp_seq.len()
                                - 1;
                        }
                    }
                }

                continue;
            }

            if rtp_loc_option.is_some() && tso_rev_comp_loc_option.is_some() {
                if rtp_loc_option.unwrap() + 1 < raw_record_len {
                    let rtp_downstream_seq =
                        if rtp_loc_option.unwrap() + tail_search_span + 1 < raw_record_len {
                            raw_record.seq()[rtp_loc_option.unwrap() + 1
                                ..=rtp_loc_option.unwrap() + tail_search_span]
                                .to_vec()
                        } else {
                            raw_record.seq()[rtp_loc_option.unwrap() + 1..].to_vec()
                        };

                    let polyt_loc_option = rtp_downstream_seq
                        .windows(tail_len)
                        .enumerate()
                        .find(|(_, window)| window.iter().all(|&b| b == 84))
                        .map(|(idx, _)| idx);
                    if polyt_loc_option.is_some() {
                        let polyt_loc = polyt_loc_option.unwrap();
                        if let Some(range) = minus_rescue_range(
                            rtp_loc_option.unwrap(),
                            polyt_loc,
                            tso_rev_comp_loc_option.unwrap(),
                            raw_record_len,
                        ) {
                            let record_number = (fq_content.len() + 1)
                                .to_string()
                                .chars()
                                .map(|c| c as u8)
                                .collect::<Vec<u8>>();

                            let rescued_record_id =
                                build_rescued_record_id(raw_record.id(), &record_number, b"minus");
                            let rescued_record_seq = raw_record.seq()[range.clone()]
                                .iter()
                                .rev()
                                .map(|c| if c & 2 == 0 { c ^ 21 } else { c ^ 4 })
                                .collect::<Vec<u8>>();
                            let rescued_record_qual = raw_record.qual().unwrap()[range]
                                .iter()
                                .rev()
                                .copied()
                                .collect::<Vec<u8>>();
                            let rescued_record = ProcessedRecord::new(
                                rescued_record_id,
                                rescued_record_seq,
                                raw_record.plus().unwrap().to_vec(),
                                rescued_record_qual,
                            );

                            fq_content.push((rescued_record, 20));
                            rescued_base_num += tso_rev_comp_loc_option.unwrap()
                                - rtp_loc_option.unwrap()
                                + tso_seq.len()
                                + rtp_seq.len()
                                - 1;
                        }
                    }
                }
            }
        }

        if !fq_content.is_empty() {
            return (fq_content, rescued_base_num);
        } else {
            let irrescuable_record_id: Vec<u8> = [b"@", raw_record.id(), b" irrescuable"].concat();
            let irrescuable_record = ProcessedRecord::new(
                irrescuable_record_id.to_owned(),
                raw_record.seq().to_vec(),
                raw_record.plus().unwrap().to_vec(),
                raw_record.qual().unwrap().to_vec(),
            );

            return (vec![(irrescuable_record, 18)], 0);
        }
    }

    let read_text_start: &[u8] = if raw_record_len > shift_threshold_start {
        &raw_record.seq()[..shift_threshold_start]
    } else {
        raw_record.seq()
    };
    let read_text_end: Vec<u8> = convert_rev(if raw_record_len > shift_threshold_end {
        &raw_record.seq()[raw_record_len - shift_threshold_end..]
    } else {
        raw_record.seq()
    });

    let (tso_ld, res_shift_len_tso) = min_edit_distance_and_coord(&tso_seq, read_text_start);
    let (rtp_ld, res_shift_len_rtp) = min_edit_distance_and_coord(&rtp_seq, read_text_start);
    let (tso_comp_ld, res_shift_len_tso_comp) =
        min_edit_distance_and_coord(&primers.tso_trimmed_comp, &read_text_end);
    let (rtp_comp_ld, res_shift_len_rtp_comp) =
        min_edit_distance_and_coord(&primers.rtp_trimmed_comp, &read_text_end);

    let res_shift_len_polya_option = read_text_end
        .windows(tail_len)
        .enumerate()
        .find(|(_, window)| window.iter().all(|&b| b == 65))
        .map(|(idx, _)| idx);
    let res_shift_len_polyt_option = read_text_start
        .windows(tail_len)
        .enumerate()
        .find(|(_, window)| window.iter().all(|&b| b == 84))
        .map(|(idx, _)| idx);

    let (res_shift_len_start, res_shift_len_end, res_fq_idx, record_id): (
        usize,
        usize,
        usize,
        Vec<u8>,
    );

    let is_high_confidence_sense_strand = tso_ld <= err_threshold_tso
        && rtp_comp_ld <= err_threshold_rtp_comp
        && res_shift_len_polya_option.is_some()
        && res_shift_len_polya_option.unwrap() >= res_shift_len_rtp_comp
        && ((umi_len == 0 && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp <= 3)
            || (umi_len > 0
                && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp
                    >= (umi_len as f64 * 0.5).round() as usize
                && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp
                    <= (umi_len as f64 * 1.5).round() as usize));
    let is_high_confidence_antisense_strand = rtp_ld <= err_threshold_rtp
        && tso_comp_ld <= err_threshold_tso_comp
        && res_shift_len_polyt_option.is_some()
        && res_shift_len_polyt_option.unwrap() >= res_shift_len_rtp
        && ((umi_len == 0 && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp <= 3)
            || (umi_len > 0
                && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp
                    >= (umi_len as f64 * 0.5).round() as usize
                && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp
                    <= (umi_len as f64 * 1.5).round() as usize));
    let is_low_confidence_sense_strand = tso_ld <= err_threshold_tso
        && res_shift_len_polya_option.is_some()
        && (rtp_comp_ld > err_threshold_rtp_comp
            || (rtp_comp_ld <= err_threshold_rtp_comp
                && (res_shift_len_polya_option.unwrap() < res_shift_len_rtp_comp
                    || (res_shift_len_polya_option.unwrap() >= res_shift_len_rtp_comp
                        && ((umi_len == 0
                            && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp
                                > 3)
                            || (umi_len > 0
                                && (res_shift_len_polya_option.unwrap()
                                    - res_shift_len_rtp_comp
                                    < (umi_len as f64 * 0.5).round() as usize
                                    || res_shift_len_polya_option.unwrap()
                                        - res_shift_len_rtp_comp
                                        > (umi_len as f64 * 1.5).round() as usize)))))));
    let is_low_confidence_antisense_strand = tso_comp_ld <= err_threshold_tso_comp
        && res_shift_len_polyt_option.is_some()
        && (rtp_ld > err_threshold_rtp
            || (rtp_ld <= err_threshold_rtp
                && (res_shift_len_polyt_option.unwrap() < res_shift_len_rtp
                    || (res_shift_len_polyt_option.unwrap() >= res_shift_len_rtp
                        && ((umi_len == 0
                            && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp > 3)
                            || (umi_len > 0
                                && (res_shift_len_polyt_option.unwrap() - res_shift_len_rtp
                                    < (umi_len as f64 * 0.5).round() as usize
                                    || res_shift_len_polyt_option.unwrap()
                                        - res_shift_len_rtp
                                        > (umi_len as f64 * 1.5).round() as usize)))))));
    let is_double_rtp_double_polyat = rtp_ld <= err_threshold_rtp
        && rtp_comp_ld <= err_threshold_rtp_comp
        && (res_shift_len_polyt_option.is_some()
            && res_shift_len_polyt_option.unwrap() >= res_shift_len_rtp
            && ((umi_len == 0 && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp <= 3)
                || (umi_len > 0
                    && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp
                        >= (umi_len as f64 * 0.5).round() as usize
                    && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp
                        <= (umi_len as f64 * 1.5).round() as usize)))
        && (res_shift_len_polya_option.is_some()
            && res_shift_len_polya_option.unwrap() >= res_shift_len_rtp_comp
            && ((umi_len == 0
                && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp <= 3)
                || (umi_len > 0
                    && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp
                        >= (umi_len as f64 * 0.5).round() as usize
                    && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp
                        <= (umi_len as f64 * 1.5).round() as usize)));
    let is_non_tso_polya_rtp = tso_ld > err_threshold_tso
        && rtp_comp_ld <= err_threshold_rtp_comp
        && (res_shift_len_polya_option.is_some()
            && res_shift_len_polya_option.unwrap() >= res_shift_len_rtp_comp
            && ((umi_len == 0
                && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp <= 3)
                || (umi_len > 0
                    && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp
                        >= (umi_len as f64 * 0.5).round() as usize
                    && res_shift_len_polya_option.unwrap() - res_shift_len_rtp_comp
                        <= (umi_len as f64 * 1.5).round() as usize)));
    let is_rtp_polyt_non_tso = rtp_ld <= err_threshold_rtp
        && tso_comp_ld > err_threshold_tso_comp
        && (res_shift_len_polyt_option.is_some()
            && res_shift_len_polyt_option.unwrap() >= res_shift_len_rtp
            && ((umi_len == 0 && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp <= 3)
                || (umi_len > 0
                    && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp
                        >= (umi_len as f64 * 0.5).round() as usize
                    && res_shift_len_polyt_option.unwrap() - res_shift_len_rtp
                        <= (umi_len as f64 * 1.5).round() as usize)));
    let is_tso_rtp_non_polyat = tso_ld <= err_threshold_tso
        && rtp_comp_ld <= err_threshold_rtp_comp
        && res_shift_len_polya_option.is_none()
        && res_shift_len_polyt_option.is_none();
    let is_rtp_tso_non_polyat = rtp_ld <= err_threshold_rtp
        && tso_comp_ld <= err_threshold_tso_comp
        && res_shift_len_polya_option.is_none()
        && res_shift_len_polyt_option.is_none();
    let is_double_tso = tso_ld <= err_threshold_tso
        && tso_comp_ld <= err_threshold_tso_comp
        && res_shift_len_polya_option.is_none()
        && res_shift_len_polyt_option.is_none();
    let is_double_rtp = rtp_ld <= err_threshold_rtp
        && rtp_comp_ld <= err_threshold_rtp_comp
        && res_shift_len_polya_option.is_none()
        && res_shift_len_polyt_option.is_none();
    let is_single_5_tso = tso_ld <= err_threshold_tso
        && rtp_comp_ld > err_threshold_rtp_comp
        && tso_comp_ld > err_threshold_tso_comp
        && rtp_ld > err_threshold_rtp;
    let is_single_5_rtp = rtp_ld <= err_threshold_rtp
        && tso_comp_ld > err_threshold_tso_comp
        && rtp_comp_ld > err_threshold_rtp_comp
        && tso_ld > err_threshold_tso;
    let is_single_3_tso = tso_comp_ld <= err_threshold_tso_comp
        && rtp_ld > err_threshold_rtp
        && tso_ld > err_threshold_tso
        && rtp_comp_ld > err_threshold_rtp_comp;
    let is_single_3_rtp = rtp_comp_ld <= err_threshold_rtp_comp
        && tso_ld > err_threshold_tso
        && rtp_ld > err_threshold_tso
        && tso_comp_ld > err_threshold_tso_comp;

    if is_high_confidence_sense_strand && !is_high_confidence_antisense_strand {
        res_shift_len_start = res_shift_len_tso;
        res_shift_len_end = if umi_len == 0 {
            res_shift_len_polya_option.unwrap()
        } else {
            res_shift_len_rtp_comp
        };
        res_fq_idx = 0;
        record_id = [b"@", raw_record.id(), b" high-confidence strand:+"].concat();
    } else if !is_high_confidence_sense_strand && is_high_confidence_antisense_strand {
        res_shift_len_start = if umi_len == 0 {
            res_shift_len_polyt_option.unwrap()
        } else {
            res_shift_len_rtp
        };
        res_shift_len_end = res_shift_len_tso_comp;
        res_fq_idx = 1;
        record_id = [b"@", raw_record.id(), b" high-confidence strand:-"].concat();
    } else if is_double_rtp_double_polyat {
        res_shift_len_start = res_shift_len_polyt_option.unwrap();
        res_shift_len_end = res_shift_len_polya_option.unwrap();
        res_fq_idx = 4;
        record_id = [b"@", raw_record.id(), b" double-rtp-double-polyat"].concat();
    } else if is_low_confidence_sense_strand {
        res_shift_len_start = res_shift_len_tso;
        res_shift_len_end = res_shift_len_polya_option.unwrap();
        res_fq_idx = 2;
        record_id = [b"@", raw_record.id(), b" low-confidence strand:+"].concat();
    } else if is_low_confidence_antisense_strand {
        res_shift_len_start = res_shift_len_polyt_option.unwrap();
        res_shift_len_end = res_shift_len_tso_comp;
        res_fq_idx = 3;
        record_id = [b"@", raw_record.id(), b" low-confidence strand:-"].concat();
    } else if is_non_tso_polya_rtp {
        res_shift_len_start = 0;
        res_shift_len_end = res_shift_len_polya_option.unwrap();
        res_fq_idx = 5;
        record_id = [b"@", raw_record.id(), b" non-tso-polya-rtp"].concat();
    } else if is_rtp_polyt_non_tso {
        res_shift_len_start = res_shift_len_polyt_option.unwrap();
        res_shift_len_end = 0;
        res_fq_idx = 6;
        record_id = [b"@", raw_record.id(), b" rtp-polyt-non-tso"].concat();
    } else if is_tso_rtp_non_polyat {
        res_shift_len_start = res_shift_len_tso;
        res_shift_len_end = res_shift_len_rtp_comp;
        res_fq_idx = 7;
        record_id = [b"@", raw_record.id(), b" tso-rtp-non-polyat"].concat();
    } else if is_rtp_tso_non_polyat {
        res_shift_len_start = res_shift_len_rtp;
        res_shift_len_end = res_shift_len_tso_comp;
        res_fq_idx = 8;
        record_id = [b"@", raw_record.id(), b" rtp-tso-non-polyat"].concat();
    } else if is_double_tso {
        res_shift_len_start = res_shift_len_tso;
        res_shift_len_end = res_shift_len_tso_comp;
        res_fq_idx = 13;
        record_id = [b"@", raw_record.id(), b" double-tso"].concat();
    } else if is_double_rtp {
        res_shift_len_start = res_shift_len_rtp;
        res_shift_len_end = res_shift_len_rtp_comp;
        res_fq_idx = 14;
        record_id = [b"@", raw_record.id(), b" double-rtp"].concat();
    } else if is_single_5_tso {
        res_shift_len_start = res_shift_len_tso;
        res_shift_len_end = 0;
        res_fq_idx = 9;
        record_id = [b"@", raw_record.id(), b" single-5'-tso"].concat();
    } else if is_single_5_rtp {
        res_shift_len_start = res_shift_len_rtp;
        res_shift_len_end = 0;
        res_fq_idx = 10;
        record_id = [b"@", raw_record.id(), b" single-5'-rtp"].concat();
    } else if is_single_3_tso {
        res_shift_len_start = 0;
        res_shift_len_end = res_shift_len_tso_comp;
        res_fq_idx = 11;
        record_id = [b"@", raw_record.id(), b" single-3'-tso"].concat();
    } else if is_single_3_rtp {
        res_shift_len_start = 0;
        res_shift_len_end = res_shift_len_rtp_comp;
        res_fq_idx = 12;
        record_id = [b"@", raw_record.id(), b" single-3'-rtp"].concat();
    } else {
        res_shift_len_start = 0;
        res_shift_len_end = 0;
        res_fq_idx = 15;
        record_id = [b"@", raw_record.id(), b" non-primer"].concat();
    }

    if raw_record_len > res_shift_len_start + res_shift_len_end {
        if res_fq_idx == 1 || res_fq_idx == 3 {
            let record_seq = raw_record.seq()
                [res_shift_len_start..raw_record_len - res_shift_len_end]
                .iter()
                .rev()
                .map(|c| if c & 2 == 0 { c ^ 21 } else { c ^ 4 })
                .collect::<Vec<u8>>();
            let record_qual = raw_record.qual().unwrap()
                [res_shift_len_start..raw_record_len - res_shift_len_end]
                .iter()
                .rev()
                .map(|&q| q)
                .collect::<Vec<u8>>();
            let splitted_record = ProcessedRecord::new(
                record_id,
                record_seq,
                raw_record.plus().unwrap().to_vec(),
                record_qual,
            );

            return (vec![(splitted_record, res_fq_idx)], raw_record_len);
        } else {
            let splitted_record = ProcessedRecord::new(
                record_id,
                raw_record.seq()[res_shift_len_start..raw_record_len - res_shift_len_end].to_vec(),
                raw_record.plus().unwrap().to_vec(),
                raw_record.qual().unwrap()[res_shift_len_start..raw_record_len - res_shift_len_end]
                    .to_vec(),
            );

            return (vec![(splitted_record, res_fq_idx)], raw_record_len);
        }
    }

    (
        vec![(
            ProcessedRecord::new(
                [b"@", raw_record.id(), b" non-primer"].concat(),
                raw_record.seq().to_vec(),
                raw_record.plus().unwrap().to_vec(),
                raw_record.qual().unwrap().to_vec(),
            ),
            15,
        )],
        0,
    )
}

#[cfg(test)]
mod tests {
    use super::{
        build_rescued_record_id, chimeric_tail_search_span, minus_rescue_range, plus_rescue_range,
        PrimerSequences,
    };
    use crate::glycine::utils::{convert_comp, convert_rev_comp};

    #[test]
    fn defaults_to_legacy_window_without_umi_gap() {
        assert_eq!(chimeric_tail_search_span(10, 0), 30);
    }

    #[test]
    fn expands_window_for_long_umi_gap() {
        assert_eq!(chimeric_tail_search_span(10, 57), 77);
    }

    #[test]
    fn rescued_id_uses_requested_underscore_format() {
        assert_eq!(
            build_rescued_record_id(b"260F200339011_1_4460_2_3736_19402_1_13.15", b"2", b"minus"),
            b"@260F200339011_1_4460_2_3736_19402_1_13.15_rescued_No2_strand_minus".to_vec()
        );
    }

    #[test]
    fn primer_variants_are_precomputed_with_legacy_transformations() {
        let tso = b"AAGACCGC".to_vec();
        let rtp = b"GAGGTCCA".to_vec();
        let trim_len = 2;
        let primers = PrimerSequences::new(tso.clone(), rtp.clone(), trim_len);

        assert_eq!(primers.tso, tso);
        assert_eq!(primers.rtp, rtp);
        assert_eq!(primers.tso_rev_comp, convert_rev_comp(&primers.tso));
        assert_eq!(primers.rtp_rev_comp, convert_rev_comp(&primers.rtp));
        assert_eq!(
            primers.tso_trimmed_comp,
            convert_comp(&primers.tso[trim_len..])
        );
        assert_eq!(
            primers.rtp_trimmed_comp,
            convert_comp(&primers.rtp[trim_len..])
        );
    }

    #[test]
    fn plus_rescue_range_rejects_inverted_or_out_of_bounds_windows() {
        assert_eq!(plus_rescue_range(10, 80, 5, 100), Some(11..75));
        assert_eq!(plus_rescue_range(79, 80, 0, 100), None);
        assert_eq!(plus_rescue_range(80, 80, 5, 100), None);
        assert_eq!(plus_rescue_range(10, 4, 5, 100), None);
        assert_eq!(plus_rescue_range(10, 120, 0, 100), None);
        assert_eq!(plus_rescue_range(usize::MAX, 80, 0, 100), None);
    }

    #[test]
    fn minus_rescue_range_rejects_inverted_or_out_of_bounds_windows() {
        assert_eq!(minus_rescue_range(10, 5, 80, 100), Some(16..80));
        assert_eq!(minus_rescue_range(79, 0, 80, 100), None);
        assert_eq!(minus_rescue_range(80, 5, 80, 100), None);
        assert_eq!(minus_rescue_range(10, 5, 120, 100), None);
        assert_eq!(minus_rescue_range(usize::MAX, 1, 80, 100), None);
    }
}
