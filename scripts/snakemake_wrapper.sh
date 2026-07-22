#!/usr/bin/env bash
set -euo pipefail

treatment=""    # Treatment sample name
control=""      # Control sample name

filter=""        # True or False for fragment length filtering
spkin=""         # True or False for spike-in calibration

threads=1

# Main output directory
out_dir=""

# Directory containing the raw FASTQ files
raw_dir=""

# FASTQ files must be gzip-compressed and named exactly as:
#   ${treatment}_R1.fastq.gz
#   ${treatment}_R2.fastq.gz
#   ${control}_R1.fastq.gz
#   ${control}_R2.fastq.gz

# The Snakefile, configuration file, helper_functions.py, and SEACR files
# are expected to be located in this directory.
pipeline_dir=""

required_vars=(
    treatment
    control
    filter
    spkin
    out_dir
    raw_dir
    pipeline_dir
)

for var in "${required_vars[@]}"; do
    if [[ -z "${!var}" ]]; then
        echo "Error: '$var' has not been set." >&2
        exit 1
    fi
done

config="$pipeline_dir/config.json"
snakefile="$pipeline_dir/Snakefile"     # or Snakefile_traditional_spkin

# Create output sub-directories
mkdir -p \
    "$out_dir/temp" \
    "$out_dir/log" \
    "$out_dir/qc" \
    "$out_dir/fastqc/pre_fastp" \
    "$out_dir/fastqc/post_fastp" \
    "$out_dir/multiqc" \
    "$out_dir/fastp" \
    "$out_dir/bam" \
    "$out_dir/bed" \
    "$out_dir/bdg" \
    "$out_dir/seacr" \
    "$out_dir/macs"

# Download SEACR 1.3 only when files are absent
SEACR_URL="https://raw.githubusercontent.com/FredHutch/SEACR/v1.3"

for file in SEACR_1.3.R SEACR_1.3.sh; do
    if [[ ! -f "$pipeline_dir/$file" ]]; then
        wget -O "$pipeline_dir/$file" "$SEACR_URL/$file"

        if [[ "$file" == "SEACR_1.3.sh" ]]; then
            chmod +x "$pipeline_dir/$file"
        fi
    fi
done

snakemake \
    --snakefile "$snakefile" \
    --cores "$threads" \
    --printshellcmds \
    -n \
    --configfile "$config" \
    --config \
        treatment="$treatment" \
        control="$control" \
        threads="$threads" \
        filter="$filter" \
        spkin="$spkin" \
        pipeline_dir="$pipeline_dir" \
        output_dir="$out_dir" \
        raw_dir="$raw_dir" \
        temp_dir="$out_dir/temp" \
        log_dir="$out_dir/log" \
        qc_dir="$out_dir/qc" \
        fastqc_dir="$out_dir/fastqc" \
        multiqc_dir="$out_dir/multiqc" \
        fastp_out_dir="$out_dir/fastp" \
        bam_dir="$out_dir/bam" \
        bed_dir="$out_dir/bed" \
        bdg_dir="$out_dir/bdg" \
        macs_dir="$out_dir/macs" \
        seacr_dir="$out_dir/seacr" \
    > "$out_dir/output.txt"
