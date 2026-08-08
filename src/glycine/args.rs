use clap::{Arg, ArgMatches, Command};
use std::ffi::OsString;
use std::path::PathBuf;

// 使用clap传参
pub fn command() -> Command {
    Command::new("glycine")
        .author("Xia Xiaoshuang <xiaxiaoshuang@genomics.cn>")
        .version("1.0.1")
        .about("A tool for identifying and trimming full-length cDNA sequencing reads")
        .args([
            Arg::new("fastq file")
                .short('f')
                .long("fastq")
                .help("Input a fastq file")
                .value_parser(clap::value_parser!(String))
                .required(true),
            Arg::new("tso seq")
                .short('5')
                .long("tso_seq")
                .help("Template switching oligo (TSO) sequence")
                .value_parser(clap::value_parser!(String))
                .required(true),
            Arg::new("rtp seq")
                .short('3')
                .long("rtp_seq")
                .help("Reverse transcription primer (RTP) sequence")
                .value_parser(clap::value_parser!(String))
                .required(true),
            Arg::new("output dir")
                .short('o')
                .long("outdir")
                .help("Output directory")
                .value_parser(clap::value_parser!(PathBuf))
                .required(true),
            Arg::new("sample")
                .short('n')
                .long("sample")
                .help("Sample name")
                .value_parser(clap::value_parser!(String))
                .required(true),
            Arg::new("err threshold")
                .short('e')
                .long("err")
                .help(
                    "Threshold of Levenshtein Distance value or sequencing error rate. \
                    Use a comma to separate two numbers, \
                    when there are different thresholds for the 5' and 3' ends",
                )
                .value_parser(clap::value_parser!(String))
                .default_value("0.25,0.25")
                .required(false),
            Arg::new("shift threshold")
                .short('s')
                .long("shift")
                .help(
                    "Threshold of shift length for identifying TSO/RTP sequences. \
                    Use a comma to separate two numbers, \
                    when there are different thresholds for the 5' and 3' ends",
                )
                .value_parser(clap::value_parser!(String))
                .default_value("100,100")
                .required(false),
            Arg::new("min len")
                .short('L')
                .long("min_len")
                .help("Sequences shorter than the minimum length will be directly classified as discarded")
                .value_parser(clap::value_parser!(usize))
                .default_value("100")
                .required(false),
            Arg::new("min qual")
                .short('Q')
                .long("min_qual")
                .help("Read quality lower than the minimum quality will be directly classified as discarded")
                .value_parser(clap::value_parser!(f64))
                .default_value("7.0")
                .required(false),
            Arg::new("trim len")
                .short('u')
                .long("trim_len")
                .help("The length of the 3' end primer sequence to be trimmed")
                .value_parser(clap::value_parser!(usize))
                .default_value("0")
                .required(false),
            Arg::new("tail len")
                .short('l')
                .long("tail_len")
                .help("PolyA tail length")
                .value_parser(clap::value_parser!(usize))
                .default_value("10")
                .required(false),
            Arg::new("umi len")
                .short('q')
                .long("umi_len")
                .help(
                    "The length of the unique molecular identifier (UMI) located between \
                    the polyA tail and the RTP sequence",
                )
                .value_parser(clap::value_parser!(usize))
                .default_value("0")
                .required(false),
            Arg::new("thread")
                .short('t')
                .long("thread")
                .help("Number of threads")
                .value_parser(clap::value_parser!(usize))
                .default_value("4")
                .required(false),
            Arg::new("keep all outputs")
                .long("keep-all-outputs")
                .help("Keep all intermediate classification FASTQ outputs instead of only the merged full-length-plus-rescued file")
                .action(clap::ArgAction::SetTrue)
                .required(false),
        ])
}

pub fn parse_argument_from<I, T>(args: I) -> ArgMatches
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    command().get_matches_from(args)
}
