use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use serde::Serialize;

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct FastqRecord {
    pub id: String,
    pub seq: String,
    pub qual: String,
}

fn reader_for(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).with_context(|| format!("open FASTQ {}", path.display()))?;
    if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

pub fn read_fastq(path: &Path) -> Result<Vec<FastqRecord>> {
    let reader = BufReader::new(reader_for(path)?);
    parse_fastq(reader)
}

pub fn for_each_fastq_batch<F>(paths: &[impl AsRef<Path>], batch_size: usize, mut callback: F) -> Result<usize>
where
    F: FnMut(Vec<FastqRecord>) -> Result<()>,
{
    let mut total = 0usize;
    let batch_size = batch_size.max(1);
    for path in paths {
        let reader = BufReader::new(reader_for(path.as_ref())?);
        total += parse_fastq_batches(reader, batch_size, &mut callback)?;
    }
    Ok(total)
}

pub fn parse_fastq<R: BufRead>(mut reader: R) -> Result<Vec<FastqRecord>> {
    let mut records = Vec::new();
    loop {
        if let Some(record) = read_one_record(&mut reader)? {
            records.push(record);
        } else {
            break;
        }
    }
    Ok(records)
}

pub fn parse_fastq_batches<R, F>(mut reader: R, batch_size: usize, callback: &mut F) -> Result<usize>
where
    R: BufRead,
    F: FnMut(Vec<FastqRecord>) -> Result<()>,
{
    let batch_size = batch_size.max(1);
    let mut total = 0usize;
    let mut batch = Vec::with_capacity(batch_size);
    while let Some(record) = read_one_record(&mut reader)? {
        batch.push(record);
        total += 1;
        if batch.len() == batch_size {
            callback(std::mem::take(&mut batch))?;
            batch = Vec::with_capacity(batch_size);
        }
    }
    if !batch.is_empty() {
        callback(batch)?;
    }
    Ok(total)
}

fn read_one_record<R: BufRead>(reader: &mut R) -> Result<Option<FastqRecord>> {
    let mut id = String::new();
    if reader.read_line(&mut id)? == 0 {
        return Ok(None);
    }
    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();
    reader.read_line(&mut seq)?;
    reader.read_line(&mut plus)?;
    reader.read_line(&mut qual)?;
    if !id.starts_with('@') || !plus.starts_with('+') {
        bail!("invalid FASTQ record near {}", id.trim());
    }
    let clean_id = id[1..].split_whitespace().next().unwrap_or("").to_string();
    Ok(Some(FastqRecord {
        id: clean_id,
        seq: seq.trim_end().to_string(),
        qual: qual.trim_end().to_string(),
    }))
}

pub fn write_fastq<W: io::Write>(mut writer: W, records: &[FastqRecord]) -> Result<()> {
    for rec in records {
        write_fastq_record(&mut writer, rec)?;
    }
    Ok(())
}

pub fn write_fastq_record<W: io::Write>(mut writer: W, rec: &FastqRecord) -> Result<()> {
    writeln!(writer, "@{}\n{}\n+\n{}", rec.id, rec.seq, rec.qual)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_basic_fastq() {
        let data = b"@r1 comment\nACGT\n+\nIIII\n@r2\nTGCA\n+\n####\n";
        let records = parse_fastq(&data[..]).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "r1");
        assert_eq!(records[1].seq, "TGCA");
    }

    #[test]
    fn streams_fastq_in_batches() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\n####\n@r3\nAAAA\n+\n!!!!\n";
        let mut batch_sizes = Vec::new();
        let total = parse_fastq_batches(&data[..], 2, &mut |batch| {
            batch_sizes.push(batch.len());
            Ok(())
        })
        .unwrap();
        assert_eq!(total, 3);
        assert_eq!(batch_sizes, vec![2, 1]);
    }
}
