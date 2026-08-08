use super::utils::ProcessedRecord;
use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::Write;

// 将拆分的record结果写入fq.gz文件
pub fn write_fq_file(
    encoder: &mut GzEncoder<File>,
    processed_record: ProcessedRecord,
) -> Result<()> {
    encoder
        .write_all(&processed_record.id)
        .with_context(|| format!("Failed to write the record.id to the file"))?;
    encoder.write_all(b"\n")?;
    encoder
        .write_all(&processed_record.seq)
        .with_context(|| format!("Failed to write the record.seq to the file"))?;
    encoder.write_all(b"\n")?;
    encoder
        .write_all(&processed_record.plus)
        .with_context(|| format!("Failed to write the record.plus to the file"))?;
    encoder.write_all(b"\n")?;
    encoder
        .write_all(&processed_record.qual)
        .with_context(|| format!("Failed to write the record.qual to the file"))?;
    encoder.write_all(b"\n")?;

    anyhow::Ok(())
}

// 将缓冲区中的数据写入底层的输出流中
pub fn flush_fq_file(res_file_vec: &mut Vec<GzEncoder<File>>) -> Result<()> {
    for encoder in res_file_vec {
        encoder
            .flush()
            .with_context(|| format!("Failed to flush the record to the file"))?;
    }

    anyhow::Ok(())
}
