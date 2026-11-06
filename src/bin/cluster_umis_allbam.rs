use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use flora::bam_runtime::bounded_hts_threads;
use flora::umi_cluster::cluster_directional;
use rayon::prelude::*;
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};
use rustc_hash::FxHashMap as HashMap;

#[derive(Debug, Parser)]
#[command(version, about = "Cluster UMIs per gene+cell and write UB-tagged BAM")]
struct Cli {
    bam: PathBuf,

    #[arg(long = "output", default_value = "tagged.sorted.bam")]
    output: PathBuf,

    #[arg(short = 'i', long = "ref_interval", default_value_t = 1000)]
    ref_interval: u32,

    #[arg(long = "cell_gene_max_reads", default_value_t = 20_000)]
    cell_gene_max_reads: usize,

    #[arg(short = 't', long = "threads", default_value_t = 4)]
    threads: usize,
}

#[derive(Clone)]
struct TagUpdate {
    umi: String,
    gene_override: Option<String>,
}

struct GroupState {
    gene: String,
    seen_reads: usize,
    umi_counts: HashMap<String, usize>,
}

struct ReadState {
    group_id: usize,
    umi: String,
    needs_gene_override: bool,
}

pub fn main() -> Result<()> {
    let cli = Cli::parse();
    let threads = bounded_hts_threads(cli.threads);
    let cluster_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("create UMI clustering thread pool")?;
    let mut bam = bam::IndexedReader::from_path(&cli.bam).with_context(|| "open indexed BAM")?;
    bam.set_threads(threads)?;
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)
        .with_context(|| format!("create {}", cli.output.display()))?;
    out.set_threads(threads)?;

    let header_view = bam.header().clone();
    let targets = header_view
        .target_names()
        .iter()
        .enumerate()
        .map(|(tid, name)| {
            (
                tid as u32,
                String::from_utf8_lossy(name).to_string(),
                header_view.target_len(tid as u32).unwrap_or(0),
            )
        })
        .collect::<Vec<_>>();

    let process_phase = Instant::now();
    for (tid, chrom, target_len) in targets {
        bam.fetch((tid, 0, target_len))?;
        let tag_updates = collect_and_cluster(&mut bam, &chrom, &cli, &cluster_pool)?;
        bam.fetch((tid, 0, target_len))?;
        for rec in bam.records() {
            let mut rec = rec?;
            if let Some(update) = tag_updates.get(rec.qname()) {
                rec.update_aux(b"UB", bam::record::Aux::String(&update.umi))?;
                if let Some(gene) = &update.gene_override {
                    rec.update_aux(b"GN", bam::record::Aux::String(gene))?;
                }
                out.write(&rec)?;
            }
        }
    }
    drop(out);
    eprintln!(
        "[timing] cluster_umis.chromosome_passes: {:.2}s",
        process_phase.elapsed().as_secs_f64()
    );
    let phase = Instant::now();
    bam::index::build(&cli.output, None, bam::index::Type::Bai, threads as u32)?;
    eprintln!(
        "[timing] cluster_umis.index: {:.2}s",
        phase.elapsed().as_secs_f64()
    );
    Ok(())
}

fn read_cap_reached(seen_reads: usize, limit: usize) -> bool {
    limit > 0 && seen_reads >= limit
}

fn collect_and_cluster(
    bam: &mut bam::IndexedReader,
    chrom: &str,
    cli: &Cli,
    cluster_pool: &rayon::ThreadPool,
) -> Result<HashMap<Vec<u8>, TagUpdate>> {
    let collect_phase = Instant::now();
    let mut group_ids: HashMap<(String, String), usize> = HashMap::default();
    let mut groups: Vec<GroupState> = Vec::new();
    let mut read_states: HashMap<Vec<u8>, ReadState> = HashMap::default();

    for rec in bam.records() {
        let rec = rec?;
        if rec.is_unmapped() {
            continue;
        }
        let rid = rec.qname().to_vec();
        let cb = get_string_tag(&rec, b"CB").context("missing CB tag")?;
        let ur = get_string_tag(&rec, b"UR").context("missing UR tag")?;
        let gn = get_string_tag(&rec, b"GN").context("missing GN tag")?;
        let (gene, needs_gene_override) = if gn == "NA" {
            (create_region_name(&rec, chrom, cli.ref_interval), true)
        } else {
            (gn, false)
        };

        let group_id = if let Some(group_id) = group_ids.get(&(cb.clone(), gene.clone())) {
            *group_id
        } else {
            let group_id = groups.len();
            group_ids.insert((cb, gene.clone()), group_id);
            groups.push(GroupState {
                gene,
                seen_reads: 0,
                umi_counts: HashMap::default(),
            });
            group_id
        };

        let group = &mut groups[group_id];
        if read_cap_reached(group.seen_reads, cli.cell_gene_max_reads) {
            continue;
        }
        group.seen_reads += 1;
        *group.umi_counts.entry(ur.clone()).or_insert(0) += 1;
        read_states.insert(
            rid,
            ReadState {
                group_id,
                umi: ur,
                needs_gene_override,
            },
        );
    }
    eprintln!(
        "[timing] cluster_umis.collect_groups[{chrom}]: {:.2}s",
        collect_phase.elapsed().as_secs_f64()
    );
    let largest = largest_group_indices(&groups, 8);
    let labels = group_ids
        .iter()
        .filter_map(|((cell, _), &index)| {
            (groups[index].umi_counts.len() >= 5_000 || largest.contains(&index))
                .then_some((index, cell.as_str()))
        })
        .collect::<HashMap<_, _>>();
    for (rank, index) in largest.into_iter().enumerate() {
        let group = &groups[index];
        let (minimum, maximum) = group
            .umi_counts
            .keys()
            .map(String::len)
            .fold((usize::MAX, 0), |(minimum, maximum), length| {
                (minimum.min(length), maximum.max(length))
            });
        eprintln!(
            "[umi_cluster] hotspot[{chrom}] rank={} cell={} gene={} reads={} distinct_umis={} length_min={} length_max={}",
            rank + 1, labels[&index], group.gene, group.seen_reads, group.umi_counts.len(),
            minimum, maximum
        );
    }

    let cluster_phase = Instant::now();
    let corrected_by_group = cluster_groups_in_pool(&groups, cluster_pool, chrom, &labels);
    eprintln!(
        "[timing] cluster_umis.cluster_groups[{chrom}]: {:.2}s",
        cluster_phase.elapsed().as_secs_f64()
    );

    let update_phase = Instant::now();
    let mut tag_updates = HashMap::default();
    for (rid, read_state) in read_states {
        if let Some(group_map) = corrected_by_group.get(read_state.group_id) {
            if let Some(umi_corr) = group_map.get(&read_state.umi) {
                let group = &groups[read_state.group_id];
                tag_updates.insert(
                    rid,
                    TagUpdate {
                        umi: umi_corr.clone(),
                        gene_override: read_state.needs_gene_override.then(|| group.gene.clone()),
                    },
                );
            }
        }
    }
    eprintln!(
        "[timing] cluster_umis.build_updates[{chrom}]: {:.2}s",
        update_phase.elapsed().as_secs_f64()
    );
    Ok(tag_updates)
}

#[cfg(test)]
fn cluster_groups(groups: &[GroupState], threads: usize) -> Result<Vec<HashMap<String, String>>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(bounded_hts_threads(threads))
        .build()
        .context("create UMI clustering thread pool")?;
    Ok(cluster_groups_in_pool(
        groups,
        &pool,
        "test",
        &HashMap::default(),
    ))
}

fn largest_group_indices(groups: &[GroupState], limit: usize) -> Vec<usize> {
    let mut largest = BinaryHeap::new();
    for (index, group) in groups.iter().enumerate() {
        largest.push(Reverse((
            group.umi_counts.len(),
            group.seen_reads,
            Reverse(index),
        )));
        if largest.len() > limit {
            largest.pop();
        }
    }
    let mut keys = largest
        .into_iter()
        .map(|Reverse(key)| key)
        .collect::<Vec<_>>();
    keys.sort_unstable_by(|a, b| b.cmp(a));
    keys.into_iter()
        .map(|(_, _, Reverse(index))| index)
        .collect()
}

fn cluster_groups_in_pool(
    groups: &[GroupState],
    pool: &rayon::ThreadPool,
    chrom: &str,
    labels: &HashMap<usize, &str>,
) -> Vec<HashMap<String, String>> {
    pool.install(|| {
        groups
            .par_iter()
            .enumerate()
            .map(|(index, group)| {
                let large = group.umi_counts.len() >= 5_000;
                let started = large.then(Instant::now);
                if large {
                    eprintln!(
                        "[umi_cluster] start[{chrom}] cell={} gene={} reads={} distinct_umis={}",
                        labels.get(&index).unwrap_or(&"unknown"),
                        group.gene,
                        group.seen_reads,
                        group.umi_counts.len()
                    );
                }
                let corrected = cluster_directional(&group.umi_counts, 3);
                if let Some(started) = started {
                    eprintln!(
                        "[umi_cluster] done[{chrom}] cell={} gene={} elapsed={:.2}s",
                        labels.get(&index).unwrap_or(&"unknown"),
                        group.gene,
                        started.elapsed().as_secs_f64()
                    );
                }
                corrected
            })
            .collect()
    })
}

fn create_region_name(rec: &bam::Record, chrom: &str, ref_interval: u32) -> String {
    let positions = rec
        .reference_positions()
        .map(|p| p as u32)
        .collect::<Vec<_>>();
    let Some(first) = positions.first().copied() else {
        return "NA".to_string();
    };
    let last = positions.last().copied().unwrap_or(first);
    region_name_from_positions(chrom, first, last, ref_interval)
}

fn region_name_from_positions(
    chrom: &str,
    start_pos: u32,
    end_pos: u32,
    ref_interval: u32,
) -> String {
    let midpoint = (start_pos + end_pos) / 2;
    let interval_start = (midpoint / ref_interval) * ref_interval;
    let interval_end = midpoint.div_ceil(ref_interval) * ref_interval;
    format!("{chrom}_{interval_start}_{interval_end}")
}

fn get_string_tag(rec: &bam::Record, tag: &[u8; 2]) -> Option<String> {
    match rec.aux(tag).ok()? {
        bam::record::Aux::String(v) => Some(v.to_string()),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::{cluster_groups, region_name_from_positions, GroupState};
    use rustc_hash::FxHashMap as HashMap;

    #[test]
    fn largest_groups_are_ranked_without_reordering_input() {
        let groups = vec![
            group("FIRST", &[("AAAA", 3)]),
            group("SECOND", &[("AAAA", 1), ("AAAT", 2)]),
            group("THIRD", &[("AAAA", 10), ("AAAT", 20)]),
            group("FOURTH", &[("AAAA", 10), ("AAAT", 20)]),
        ];
        assert_eq!(super::largest_group_indices(&groups, 2), vec![2, 3]);
        assert_eq!(
            super::largest_group_indices(&groups, 0),
            Vec::<usize>::new()
        );
        assert_eq!(groups[0].gene, "FIRST");
        assert_eq!(super::largest_group_indices(&[], 5), Vec::<usize>::new());
    }

    #[test]
    fn zero_read_cap_means_unlimited() {
        assert!(!super::read_cap_reached(50000, 0));
        assert!(super::read_cap_reached(20000, 20000));
        assert!(!super::read_cap_reached(19999, 20000));
    }

    #[test]
    fn region_name_matches_python_midpoint_binning() {
        assert_eq!(
            region_name_from_positions("chr1", 100, 199, 1000),
            "chr1_0_1000"
        );
        assert_eq!(
            region_name_from_positions("chr1", 1000, 1000, 1000),
            "chr1_1000_1000"
        );
        assert_eq!(
            region_name_from_positions("chr7", 1444, 1555, 1000),
            "chr7_1000_2000"
        );
    }

    #[test]
    fn parallel_group_clustering_matches_single_thread_exactly() {
        let groups = vec![
            group("GENE1", &[("AAAA", 20), ("AAAT", 3), ("AATT", 1)]),
            group("GENE2", &[("CCCC", 8), ("CCCA", 2), ("GGGG", 7)]),
            group("GENE3", &[("TTTT", 1)]),
        ];
        assert_eq!(
            cluster_groups(&groups, 1).unwrap(),
            cluster_groups(&groups, 4).unwrap()
        );
    }

    fn group(gene: &str, counts: &[(&str, usize)]) -> GroupState {
        let mut umi_counts = HashMap::default();
        for (umi, count) in counts {
            umi_counts.insert((*umi).to_string(), *count);
        }
        GroupState {
            gene: gene.to_string(),
            seen_reads: counts.iter().map(|(_, count)| count).sum(),
            umi_counts,
        }
    }
}
