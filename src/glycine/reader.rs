use anyhow::{self, Context, Result};
use fxread::{self, FastxRead, Record};

// 读取输入的fastq文件
pub fn read_fq_file(fq_file_path: &String) -> Result<Box<dyn FastxRead<Item = Record>>> {
    let fq_reader = fxread::initialize_reader(fq_file_path)
        .with_context(|| format!("Failed to open the input fastq file {}", fq_file_path))?;
    anyhow::Ok(fq_reader)
}
