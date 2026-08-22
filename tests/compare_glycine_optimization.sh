#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash tests/compare_glycine_optimization.sh \
    OLD_FLORA NEW_FLORA INPUT.fastq.gz WORK_DIR TSO_SEQ RTP_SEQ [BENCH_THREADS]

The script first runs both release binaries with one thread and compares:
  1. identifying_statistic.txt byte for byte
  2. every gzip FASTQ after decompression byte for byte

If BENCH_THREADS is supplied, it then reruns both binaries with that same
thread count and writes POSIX wall-time summaries below WORK_DIR/benchmark.
EOF
}

if (( $# < 6 || $# > 7 )); then
  usage >&2
  exit 2
fi

old_flora=$1
new_flora=$2
input_fastq=$3
work_dir=$4
tso_seq=$5
rtp_seq=$6
bench_threads=${7:-}

for path in "$old_flora" "$new_flora" "$input_fastq"; do
  if [[ ! -e "$path" ]]; then
    echo "Missing input: $path" >&2
    exit 1
  fi
done

old_flora=$(cd "$(dirname "$old_flora")" && pwd)/$(basename "$old_flora")
new_flora=$(cd "$(dirname "$new_flora")" && pwd)/$(basename "$new_flora")
input_fastq=$(cd "$(dirname "$input_fastq")" && pwd)/$(basename "$input_fastq")
mkdir -p "$work_dir"
work_dir=$(cd "$work_dir" && pwd)

run_glycine() {
  local binary=$1
  local output=$2
  local threads=$3
  rm -rf "$output"
  mkdir -p "$output"
  "$binary" glycine \
    --fastq "$input_fastq" \
    --tso_seq "$tso_seq" \
    --rtp_seq "$rtp_seq" \
    --outdir "$output" \
    --sample parity \
    --jobs 1 \
    --total-threads "$threads" \
    --keep-all-outputs
}

old_correctness="$work_dir/correctness/old"
new_correctness="$work_dir/correctness/new"
run_glycine "$old_flora" "$old_correctness" 1
run_glycine "$new_flora" "$new_correctness" 1

cmp \
  "$old_correctness/parity.identifying_statistic.txt" \
  "$new_correctness/parity.identifying_statistic.txt"

old_fastq_list="$work_dir/correctness/old-fastq-files.txt"
new_fastq_list="$work_dir/correctness/new-fastq-files.txt"
(
  cd "$old_correctness"
  find . -maxdepth 1 -type f \( -name '*.fq.gz' -o -name '*.fastq.gz' \) -print |
    sed 's#^\./##' | sort
) >"$old_fastq_list"
(
  cd "$new_correctness"
  find . -maxdepth 1 -type f \( -name '*.fq.gz' -o -name '*.fastq.gz' \) -print |
    sed 's#^\./##' | sort
) >"$new_fastq_list"

if ! cmp "$old_fastq_list" "$new_fastq_list"; then
  echo "FASTQ output file lists differ" >&2
  diff -u "$old_fastq_list" "$new_fastq_list" >&2 || true
  exit 1
fi

while IFS= read -r name; do
  cmp \
    <(gzip -dc "$old_correctness/$name") \
    <(gzip -dc "$new_correctness/$name")
done <"$old_fastq_list"

echo "PASS: single-thread statistics and decompressed FASTQ outputs are identical"

if [[ -n "$bench_threads" ]]; then
  if (( bench_threads < 1 )); then
    echo "BENCH_THREADS must be greater than zero" >&2
    exit 2
  fi
  benchmark_dir="$work_dir/benchmark"
  mkdir -p "$benchmark_dir"
  rm -rf "$benchmark_dir/old" "$benchmark_dir/new"
  mkdir -p "$benchmark_dir/old" "$benchmark_dir/new"

  /usr/bin/time -p -o "$benchmark_dir/old.time.txt" \
    "$old_flora" glycine \
      --fastq "$input_fastq" \
      --tso_seq "$tso_seq" \
      --rtp_seq "$rtp_seq" \
      --outdir "$benchmark_dir/old" \
      --sample benchmark \
      --jobs 1 \
      --total-threads "$bench_threads"
  /usr/bin/time -p -o "$benchmark_dir/new.time.txt" \
    "$new_flora" glycine \
      --fastq "$input_fastq" \
      --tso_seq "$tso_seq" \
      --rtp_seq "$rtp_seq" \
      --outdir "$benchmark_dir/new" \
      --sample benchmark \
      --jobs 1 \
      --total-threads "$bench_threads"

  echo "Old release timing:"
  cat "$benchmark_dir/old.time.txt"
  echo "New release timing:"
  cat "$benchmark_dir/new.time.txt"
fi
