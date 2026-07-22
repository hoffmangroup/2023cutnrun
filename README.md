# CUT&RUN Analysis Pipeline

This repository contains the Snakemake workflows and supporting scripts used to reproduce the analyses described in the paper:

> **"Benchmarking CUT&RUN analysis using motif enrichment"**

The repository includes two workflows:

- **`Snakefile`** – the primary benchmarking pipeline used throughout the paper.
- **`Snakefile_traditional_spkin`** – an alternative workflow implementing traditional spike-in calibration using an arbitrary scale.

---

# Required Files

To run the primary workflow (`Snakefile`), the following files must be located in the same directory:

- `Snakefile`
- `config.json`
- `helper_functions.py`
- `multiqc_config.yaml`

The following external files are also required, and their locations should be specified in `config.json`:

| File                                                      | Purpose                                                             |
| --------------------------------------------------------- | ------------------------------------------------------------------- |
| `JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt` | JASPAR motif database used for motif enrichment analysis.           |
| `spike_in_calibration_new.csh`                            | Empirical spike-in calibration script used by the primary workflow. |

To reproduce the traditional spike-in calibration analysis (`Snakefile_traditional_spkin`), the original `spike_in_calibration.csh` script is required. This script is not included in this repository but is available from the Henikoff laboratory GitHub repository:

https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.csh

---

# Running the Pipeline

The primary workflow is executed using

```bash
bash snakemake_wrapper.sh
```

Before running the pipeline, edit `snakemake_wrapper.sh` and specify

- `treatment`
- `control`
- `filter`
- `spkin`
- `threads`
- `raw_dir`
- `out_dir`

The wrapper script automatically

- creates the required output directories,
- downloads SEACR v1.3 (if not already present),
- launches the Snakemake workflow.

---

# Input

The pipeline can be started from multiple stages of the workflow. If intermediate files are already available, Snakemake automatically skips the corresponding upstream steps.

| Starting point      | Expected files                                    |
| ------------------- | ------------------------------------------------- |
| Raw sequencing data | `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz` |
| Alignment           | `<sample>_main.bam`                               |

For paired analyses, `<sample>` should be replaced with the treatment or control sample name specified in `snakemake_wrapper.sh`.

---

# Reproducing the Analyses in the Manuscript

In the manuscript, each sample was processed four times to evaluate all combinations of fragment length filtering and spike-in calibration:

| Fragment length filtering | Spike-in calibration |
| ------------------------- | -------------------- |
| No                        | No                   |
| Yes                       | No                   |
| No                        | Yes                  |
| Yes                       | Yes                  |

To avoid repeating the computationally intensive preprocessing steps, we recommend the following workflow:

1. Run the pipeline once starting from the raw paired-end FASTQ files.
2. For each subsequent analysis, create symbolic links to the previously generated:
   - BAM files (`bam_dir`)
   - FastQC results (`fastqc_dir`)
   - fastp output files (`fastp_out_dir`)
3. Start the pipeline from the BAM files (`<sample>_main.bam`) for the remaining analyses.

This workflow reproduces the analyses presented in the manuscript while avoiding redundant alignment and preprocessing.

---

# Configuration File

The pipeline is configured through `config.json`.

When the pipeline is run using `snakemake_wrapper.sh`, the wrapper automatically sets the sample names and directory paths.

Before running the pipeline, the following fields must be filled in:

- `bowtie2_index_main`
- `bowtie2_index_spike_in`
- `genome_main`
- `genome_size`
- `chr_len_file`

Other parameters should also be reviewed and changed as appropriate for the samples.

## Sample Information

| Parameter   | Description            |
| ----------- | ---------------------- |
| `treatment` | Treatment sample name. |
| `control`   | Control sample name.   |

## Input and Output

| Parameter       | Description                                   |
| --------------- | --------------------------------------------- |
| `raw_dir`       | Directory containing the raw FASTQ files.     |
| `output_dir`    | Main output directory.                        |
| `temp_dir`      | Directory containing temporary files.         |
| `log_dir`       | Directory containing log files.               |
| `qc_dir`        | Directory containing quality control results. |
| `fastqc_dir`    | Directory containing FastQC results.          |
| `multiqc_dir`   | Directory containing MultiQC reports.         |
| `fastp_out_dir` | Directory containing fastp output files.      |
| `bam_dir`       | Directory containing BAM files.               |
| `bed_dir`       | Directory containing BED files.               |
| `bdg_dir`       | Directory containing bedGraph files.          |
| `macs_dir`      | Directory containing MACS output files.       |
| `seacr_dir`     | Directory containing SEACR output files.      |

## General Parameter

| Parameter | Description            |
| --------- | ---------------------- |
| `threads` | Number of CPU threads. |

## Preprocessing Parameter

| Parameter     | Description                   |
| ------------- | ----------------------------- |
| `fastp_param` | Parameters passed to `fastp`. |

## Alignment

| Parameter                | Description                                                                  |
| ------------------------ | ---------------------------------------------------------------------------- |
| `bowtie2_main_param`     | Parameters passed to Bowtie2 for alignment to the target reference genome.   |
| `bowtie2_spike_in_param` | Parameters passed to Bowtie2 for alignment to the spike-in reference genome. |
| `bowtie2_index_main`     | Path and prefix of the Bowtie2 index for the target reference genome.        |
| `bowtie2_index_spike_in` | Path and prefix of the Bowtie2 index for the spike-in reference genome.      |

**Note:** The Bowtie2 index should be specified as the path and filename prefix, **not** the directory or individual index files. For example, if the index files are

```
hg38/Bowtie2Index/
├── genome.1.bt2
├── genome.2.bt2
├── genome.3.bt2
├── genome.4.bt2
├── genome.rev.1.bt2
└── genome.rev.2.bt2
```

then `bowtie2_index_main` should be set to

```
hg38/Bowtie2Index/genome
```

## Reference Genome

| Parameter           | Description                                                                  |
| ------------------- | ---------------------------------------------------------------------------- |
| `genome_main`       | Path to the target reference genome FASTA file.                              |
| `genome_size`       | Genome size of the target reference genome.                                  |
| `chr_len_file`      | Path to the chromosome-length (`.fai`) file for the target reference genome. |
| `qualimap_assembly` | Genome assembly used by Qualimap.                                            |

## Fragment Length Filtering

| Parameter       | Description                               |
| --------------- | ----------------------------------------- |
| `filter`        | Enable fragment length filtering.         |
| `max_frag_len`  | Maximum fragment length.                  |
| `min_frag_len`  | Minimum fragment length.                  |
| `filter_option` | Fragment filtering mode (`min` or `max`). |

## Spike-in Calibration

| Parameter          | Description                                                                      |
| ------------------ | -------------------------------------------------------------------------------- |
| `spkin`            | Enable spike-in calibration.                                                     |
| `spike_in_script`  | Path to the SEACR spike-in calibration script.                                   |
| `spike_in_min_len` | Minimum fragment length used to generate the SEACR spike-in-normalized bedGraph. |
| `spike_in_max_len` | Maximum fragment length used to generate the SEACR spike-in-normalized bedGraph. |
| `scale`            | Arbitrary scale used by the traditional spike-in calibration workflow.           |

## Peak Calling

| Parameter     | Description                                     |
| ------------- | ----------------------------------------------- |
| `peak_caller` | Peak caller to use (`MACS`, `SEACR`, or `all`). |

### MACS

| Parameter         | Description                                                                              |
| ----------------- | ---------------------------------------------------------------------------------------- |
| `macs_version`    | MACS executable.                                                                         |
| `macs_fdr`        | False discovery rate threshold used by MACS.                                             |
| `min_peak_length` | Minimum retained MACS peak length.                                                       |
| `macs_callpeak`   | Use the `macs2 callpeak` command instead of the MACS subcommands used in the manuscript. |
| `cutoff`          | Perform MACS cutoff analysis.                                                            |

### SEACR

| Parameter      | Description                            |
| -------------- | -------------------------------------- |
| `seacr_script` | Path to the SEACR executable.          |
| `seacr_mode`   | SEACR mode (`stringent` or `relaxed`). |

## Motif Analysis

| Parameter       | Description                                                         |
| --------------- | ------------------------------------------------------------------- |
| `motif`         | Perform motif analysis.                                             |
| `jaspar`        | Path to the JASPAR motif database in MEME format.                   |
| `peak_location` | Peak location used for motif analysis (`summit`, `peak`, or `all`). |
| `fasta_ext_len` | Extension length used when extracting sequences for motif analysis. |
| `skip_dreme`    | Skip the MEME, DREME, FIMO, and SpaMo analyses.                     |

## Quality Control and Reporting

| Parameter        | Description                                              |
| ---------------- | -------------------------------------------------------- |
| `flagstat`       | Generate alignment statistics using Sambamba `flagstat`. |
| `fastqc`         | Run FastQC.                                              |
| `qc`             | Generate quality control metrics.                        |
| `manual`         | Only generate the MultiQC report from existing results.  |
| `multiqc_config` | Path to the MultiQC configuration file.                  |
