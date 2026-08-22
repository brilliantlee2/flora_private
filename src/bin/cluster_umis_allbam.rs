use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use flora::bam_runtime::bounded_hts_threads;
use flora::umi_cluster::cluster_directional;
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
        let tag_updates = collect_and_cluster(&mut bam, &chrom, &cli)?;
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

fn collect_and_cluster(
    bam: &mut bam::IndexedReader,
    chrom: &str,
    cli: &Cli,
) -> Result<HashMap<Vec<u8>, TagUpdate>> {
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
        if group.seen_reads >= cli.cell_gene_max_reads {
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

    let corrected_by_group = groups
        .iter()
        .map(|group| cluster_directional(&group.umi_counts, 3))
        .collect::<Vec<_>>();

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
    Ok(tag_updates)
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
    use super::region_name_from_positions;

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
}
