#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./split_and_run.sh <input_file> <chunk_size> [PREFIX] [LOG_FILE]
#
# Example:
#   ./split_and_run.sh biglist.txt 10 chunk LOG.log

input="${1:-}"
chunk_size="${2:-}"
prefix="${3:chunk}"
LOG_FILE="${4:-LOG.log}"

# Optional: keep chunk/logs in a dedicated folder
chunk_dir="."
log_dir="."

if [[ -z "$input" || -z "$chunk_size" ]]; then
  echo "Usage: $0 <input_file> <chunk_size> [PREFIX] [LOG_FILE]"
  exit 1
fi
if [[ ! -f "$input" ]]; then
  echo "Error: input file not found: $input"
  exit 1
fi
if ! [[ "$chunk_size" =~ ^[0-9]+$ ]] || [[ "$chunk_size" -le 0 ]]; then
  echo "Error: chunk_size must be a positive integer (got: $chunk_size)"
  exit 1
fi

total_lines=$(wc -l < "$input")
n_chunks=$(( (total_lines + chunk_size - 1) / chunk_size ))

echo "=== $(date '+%F %T') : starting (input=$input, chunk_size=$chunk_size, prefix=$prefix) ===" | tee -a "$LOG_FILE"

echo "Input file      : $input"        | tee -a "$LOG_FILE"
echo "Total lines     : $total_lines"  | tee -a "$LOG_FILE"
echo "Chunk size      : $chunk_size"   | tee -a "$LOG_FILE"
echo "Expected chunks : $n_chunks"     | tee -a "$LOG_FILE"


count=1
while true; do
  start=$(( (count - 1) * chunk_size + 1 ))
  end=$(( start + chunk_size - 1 ))

  chunk_file="${chunk_dir}/${prefix}_${count}.txt"
  stdout_log="${log_dir}/root_log_${count}.log"
  stderr_log="${log_dir}/root_log_${count}.err"

  # Create chunk file with output file line + the slice of the input
  echo "output_${prefix}_${count}.root" > "$chunk_file"
  sed -n "${start},${end}p" "$input" >> "$chunk_file"

  # If the extracted slice is empty -> chunk file will only contain the output file line
  if [[ $(wc -l < "$chunk_file") -eq 1 ]]; then
    rm -f "$chunk_file"
    break
  fi

  echo "$(date '+%F %T') [RUN ] chunk $count (lines ${start}-${end}) using $chunk_file" | tee -a "$LOG_FILE"

  if root -b -q "epic_studies.cpp(\"$chunk_file\")" >"$stdout_log" 2>"$stderr_log"; then
    echo "$(date '+%F %T') [ OK ] chunk $count" | tee -a "$LOG_FILE"
  else
    rc=$?
    echo "$(date '+%F %T') [FAIL] chunk $count (exit code $rc)" | tee -a "$LOG_FILE"
    # If you want to stop at first failure, uncomment:
    # exit "$rc"
  fi

  count=$((count + 1))
done

echo "=== $(date '+%F %T') : all done (processed $((count-1)) chunks) ===" | tee -a "$LOG_FILE"
