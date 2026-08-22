use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use csv::{ReaderBuilder, StringRecord, WriterBuilder};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

const FRACTIONS: [f64; 13] = [
    0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
];

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Compute sequencing saturation from a cell/read table"
)]
struct Cli {
    #[arg(long, default_value = "cell_umi_gene.tsv")]
    input: PathBuf,
    #[arg(long = "output-tsv", default_value = "saturation.tsv")]
    output_tsv: PathBuf,
    // Accepted for workflow CLI compatibility; plotting remains in Python.
    #[arg(long = "output-png")]
    _output_png: Option<PathBuf>,
}

#[derive(Clone, Copy)]
struct Row {
    read: u32,
    gene: u32,
    cell: u32,
    umi: u32,
    molecule: u32,
    known_gene: bool,
}

pub fn main() -> Result<()> {
    let cli = Cli::parse();
    let (mut rows, n_cells) = load_rows(&cli.input)?;
    numpy_random_state_42_shuffle(&mut rows);
    let records = compute_records(&rows, n_cells, &FRACTIONS);
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&cli.output_tsv)?;
    writer.write_record([
        "fraction",
        "reads",
        "reads_per_cell",
        "genes_per_cell",
        "umis_per_cell",
        "saturation",
    ])?;
    for row in records {
        writer.serialize(row)?;
    }
    writer.flush()?;
    Ok(())
}

fn load_rows(path: &PathBuf) -> Result<(Vec<Row>, usize)> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("open {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let index = |name: &str| {
        headers
            .iter()
            .position(|value| value == name)
            .with_context(|| format!("missing {name} column"))
    };
    let read_idx = index("read_id")?;
    let gene_idx = index("gene")?;
    let cell_idx = index("barcode")?;
    let umi_idx = index("umi")?;
    let mut reads = HashMap::default();
    let mut genes = HashMap::default();
    let mut cells = HashMap::default();
    let mut umis = HashMap::default();
    let mut molecules = HashMap::default();
    let mut rows = Vec::new();
    for record in reader.records() {
        let record = record?;
        let read = text(&record, read_idx);
        let gene = text(&record, gene_idx);
        let cell = text(&record, cell_idx);
        let umi = text(&record, umi_idx);
        // Preserve the historical Python concatenation semantics exactly.
        let molecule = format!("{gene}_{cell}_{umi}");
        rows.push(Row {
            read: intern(&mut reads, read),
            gene: intern(&mut genes, gene),
            cell: intern(&mut cells, cell),
            umi: intern(&mut umis, umi),
            molecule: intern(&mut molecules, &molecule),
            known_gene: is_known_gene(gene),
        });
    }
    Ok((rows, cells.len()))
}

fn text(record: &StringRecord, index: usize) -> &str {
    record.get(index).unwrap_or("")
}

fn intern(map: &mut HashMap<String, u32>, value: &str) -> u32 {
    if let Some(id) = map.get(value) {
        return *id;
    }
    let id = map.len() as u32;
    map.insert(value.to_string(), id);
    id
}

fn is_known_gene(value: &str) -> bool {
    let value = value.trim();
    if matches!(value, "" | "NA" | "nan" | "None") {
        return false;
    }
    !looks_like_genomic_placeholder(value)
}

fn looks_like_genomic_placeholder(value: &str) -> bool {
    let bytes = value.as_bytes();
    for first in 0..bytes.len() {
        if !bytes[first].is_ascii_alphanumeric() {
            continue;
        }
        let mut i = first;
        while i < bytes.len() && bytes[i].is_ascii_alphanumeric() {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'_' {
            continue;
        }
        i += 1;
        let start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i == start || i >= bytes.len() || bytes[i] != b'_' {
            continue;
        }
        i += 1;
        let start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i > start {
            return true;
        }
    }
    false
}

fn compute_records(
    rows: &[Row],
    n_cells: usize,
    fractions: &[f64],
) -> Vec<(f64, usize, f64, f64, f64, f64)> {
    let sizes: Vec<usize> = fractions
        .iter()
        .map(|fraction| round_half_even(*fraction * rows.len() as f64) as usize)
        .collect();
    let mut read_seen = HashSet::default();
    let mut umi_seen = HashSet::default();
    let mut gene_seen = HashSet::default();
    let mut molecule_seen = HashSet::default();
    let mut active = vec![false; n_cells];
    let mut read_counts = vec![0usize; n_cells];
    let mut umi_counts = vec![0usize; n_cells];
    let mut gene_counts = vec![0usize; n_cells];
    let mut output = Vec::with_capacity(fractions.len());
    let mut consumed = 0usize;
    for (&fraction, &size) in fractions.iter().zip(&sizes) {
        while consumed < size {
            let row = rows[consumed];
            let cell = row.cell as usize;
            active[cell] = true;
            if read_seen.insert(pair(row.cell, row.read)) {
                read_counts[cell] += 1;
            }
            if umi_seen.insert(pair(row.cell, row.umi)) {
                umi_counts[cell] += 1;
            }
            if row.known_gene && gene_seen.insert(pair(row.cell, row.gene)) {
                gene_counts[cell] += 1;
            }
            molecule_seen.insert(row.molecule);
            consumed += 1;
        }
        if size == 0 {
            output.push((fraction, 0, 0.0, 0.0, 0.0, 0.0));
            continue;
        }
        output.push((
            fraction,
            size,
            median_active(&read_counts, &active),
            median_active(&gene_counts, &active),
            median_active(&umi_counts, &active),
            1.0 - molecule_seen.len() as f64 / size as f64,
        ));
    }
    output
}

fn pair(left: u32, right: u32) -> u64 {
    ((left as u64) << 32) | right as u64
}

fn median_active(values: &[usize], active: &[bool]) -> f64 {
    let mut selected: Vec<usize> = values
        .iter()
        .zip(active)
        .filter_map(|(value, is_active)| is_active.then_some(*value))
        .collect();
    selected.sort_unstable();
    let middle = selected.len() / 2;
    if selected.len() % 2 == 0 {
        (selected[middle - 1] + selected[middle]) as f64 / 2.0
    } else {
        selected[middle] as f64
    }
}

fn round_half_even(value: f64) -> u64 {
    let floor = value.floor();
    let fraction = value - floor;
    if fraction < 0.5 {
        floor as u64
    } else if fraction > 0.5 {
        floor as u64 + 1
    } else if floor as u64 % 2 == 0 {
        floor as u64
    } else {
        floor as u64 + 1
    }
}

struct Mt19937 {
    state: [u32; 624],
    index: usize,
}

impl Mt19937 {
    fn seeded(seed: u32) -> Self {
        let mut state = [0; 624];
        state[0] = seed;
        for i in 1..624 {
            state[i] = 1812433253u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, index: 624 }
    }

    fn next(&mut self) -> u32 {
        if self.index >= 624 {
            self.twist();
        }
        let mut value = self.state[self.index];
        self.index += 1;
        value ^= value >> 11;
        value ^= (value << 7) & 0x9d2c5680;
        value ^= (value << 15) & 0xefc60000;
        value ^= value >> 18;
        value
    }

    fn twist(&mut self) {
        for i in 0..624 {
            let value = (self.state[i] & 0x80000000) | (self.state[(i + 1) % 624] & 0x7fffffff);
            self.state[i] = self.state[(i + 397) % 624]
                ^ (value >> 1)
                ^ if value & 1 != 0 { 0x9908b0df } else { 0 };
        }
        self.index = 0;
    }

    fn interval(&mut self, max: u32) -> u32 {
        let mask = if max == 0 {
            0
        } else {
            u32::MAX >> max.leading_zeros()
        };
        loop {
            let value = self.next() & mask;
            if value <= max {
                return value;
            }
        }
    }
}

fn numpy_random_state_42_shuffle<T>(values: &mut [T]) {
    let mut rng = Mt19937::seeded(42);
    for i in (1..values.len()).rev() {
        let j = rng.interval(i as u32) as usize;
        values.swap(i, j);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numpy_seed_42_permutation_matches_reference() {
        let mut values: Vec<_> = (0..8).collect();
        numpy_random_state_42_shuffle(&mut values);
        assert_eq!(values, vec![1, 5, 0, 7, 2, 4, 3, 6]);
    }

    #[test]
    fn half_even_rounding_matches_pandas_sample_size() {
        assert_eq!(round_half_even(2.5), 2);
        assert_eq!(round_half_even(3.5), 4);
    }

    #[test]
    fn underscore_molecule_collision_is_preserved() {
        let rows = vec![
            Row {
                read: 0,
                gene: 0,
                cell: 0,
                umi: 0,
                molecule: 0,
                known_gene: true,
            },
            Row {
                read: 1,
                gene: 1,
                cell: 1,
                umi: 0,
                molecule: 0,
                known_gene: true,
            },
        ];
        let result = compute_records(&rows, 2, &[1.0]);
        assert_eq!(result[0].5, 0.5);
    }
}
