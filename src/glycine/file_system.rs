use anyhow::{self, Context, Result};
use flate2::{write::GzEncoder, Compression};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

// 创建结果输出目录和文件
pub fn create_output_dir(
    output_dir_path: &PathBuf,
    sample: &String,
    keep_all_outputs: bool,
) -> Result<Vec<GzEncoder<File>>> {
    if !output_dir_path.exists() {
        fs::create_dir_all(output_dir_path).with_context(|| {
            format!(
                "Failed to create the directory {}",
                output_dir_path.display()
            )
        })?;
    }

    let full_length_fq_path = output_dir_path.join(format!("{}.full-length.fq.gz", sample));
    let non_full_length_fq_path = output_dir_path.join(format!("{}.non-full-length.fq.gz", sample));
    let discarded_fq_path = output_dir_path.join(format!("{}.discarded.fq.gz", sample));
    let failed_filter_fq_path = output_dir_path.join(format!("{}.failed-filter.fq.gz", sample));
    let rescued_fq_path = output_dir_path.join(format!("{}.rescued.fq.gz", sample));
    let merged_fq_path = output_dir_path.join(format!("{}.full-length-plus-rescued.fq.gz", sample));

    let output_files_path = if keep_all_outputs {
        vec![
            full_length_fq_path,
            non_full_length_fq_path,
            discarded_fq_path,
            failed_filter_fq_path,
            rescued_fq_path,
            merged_fq_path,
        ]
    } else {
        vec![merged_fq_path]
    };
    let mut encoders = Vec::new();

    for file_path in output_files_path {
        let file = File::create(&file_path)
            .with_context(|| format!("Failed to create the file {}", file_path.display()))?;

        let encoder = GzEncoder::new(file, Compression::default());
        encoders.push(encoder);
    }

    anyhow::Ok(encoders)
}

// 创建identifying_statistic.txt文件并将拆分的统计结果写入
pub fn create_and_write_statistic_file(
    output_dir_path: &PathBuf,
    sample: &String,
    content: String,
) -> Result<()> {
    let statistic_file_path = output_dir_path.join(format!("{}.identifying_statistic.txt", sample));

    let mut statistic_file = File::create(&statistic_file_path).with_context(|| {
        format!(
            "Failed to create the file {}",
            statistic_file_path.display()
        )
    })?;

    statistic_file
        .write_all(content.as_bytes())
        .with_context(|| format!("Failed to write the file {}", statistic_file_path.display()))?;

    anyhow::Ok(())
}

#[cfg(test)]
mod tests {
    use super::create_output_dir;
    use std::path::PathBuf;

    #[test]
    fn create_output_dir_creates_merged_fastq_output() {
        let temp_dir =
            std::env::temp_dir().join(format!("glycine-output-test-{}", std::process::id()));
        let sample = String::from("demo");

        let encoders = create_output_dir(&PathBuf::from(&temp_dir), &sample, true).unwrap();

        assert_eq!(encoders.len(), 6);
        assert!(temp_dir
            .join("demo.full-length-plus-rescued.fq.gz")
            .exists());

        drop(encoders);
        std::fs::remove_dir_all(temp_dir).unwrap();
    }

    #[test]
    fn create_output_dir_compact_mode_only_creates_merged_output() {
        let temp_dir = std::env::temp_dir().join(format!(
            "glycine-output-compact-test-{}",
            std::process::id()
        ));
        let sample = String::from("demo");

        let encoders = create_output_dir(&PathBuf::from(&temp_dir), &sample, false).unwrap();

        assert_eq!(encoders.len(), 1);
        assert!(temp_dir
            .join("demo.full-length-plus-rescued.fq.gz")
            .exists());
        assert!(!temp_dir.join("demo.full-length.fq.gz").exists());

        drop(encoders);
        std::fs::remove_dir_all(temp_dir).unwrap();
    }
}
