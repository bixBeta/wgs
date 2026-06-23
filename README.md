[![Build and Push Docker Images](https://github.com/bixBeta/wgs/actions/workflows/docker-build.yml/badge.svg)](https://github.com/bixBeta/wgs/actions/workflows/docker-build.yml)

# WGS Workflow — BioHPC Cornell

Nextflow DSL2 pipeline for whole-genome sequencing analysis on BioHPC Cornell servers.

**Modules:** fastp → bowtie2 → samtools → picard (MarkDuplicates) → qualimap → deeptools (GC bias) → MultiQC → Quarto GC bias report

---

## Prerequisites (BioHPC)

- **Nextflow** ≥ 23.x (`module load nextflow` or available in your PATH)
- **Apptainer/Singularity** — available by default on BioHPC nodes (Docker is not required)
- No local software installation required — all tools run inside Singularity containers pulled automatically from GHCR

---

## Quick Start

```bash
nextflow run bixBeta/wgs -r container \
  --id PROJECT_ID \
  --sheet sample-sheet.csv \
  --genome hg38 \
  --bowtie2 \
  --fastp
```

---

## Sample Sheet

Tab or comma-separated file with three columns:

```
label,fastq1,fastq2
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz
SAMPLE2,/path/to/SAMPLE2_R1.fastq.gz,/path/to/SAMPLE2_R2.fastq.gz
```

> If your FASTQs are in a `fastqs/` subdirectory in the run folder, pass `--fastqs` and use just the filenames (no path) in the sheet.

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--id` | `TREx_ID` | Project ID — used for output file naming |
| `--sheet` | `sample-sheet.csv` | Path to sample sheet |
| `--genome` | `null` | Reference genome (see below) |
| `--mode` | `PE` | Sequencing mode: `PE` (paired-end) |
| `--bowtie2` | `false` | Enable bowtie2 alignment |
| `--fastp` | `false` | Enable fastp trimming (recommended) |
| `--gcbias` | `false` | Enable GC bias analysis (requires `--bowtie2`) |
| `--fastqs` | `false` | Look for FASTQs in `./fastqs/` subdirectory |
| `--listGenomes` | `false` | Print available genome index paths and exit |

---

## Available Genomes

Bowtie2 indices live at `/workdir/genomes/` on the BioHPC server.

| Key | Species | Assembly |
|-----|---------|----------|
| `hg38` | Human | GRCh38 |
| `mm10` | Mouse | GRCm38 |
| `dm6` | Drosophila | BDGP6 |
| `canFam4` | Dog | canFam4 |
| `fc9` | Cat | Felis_catus_9.0 |
| `vitis` | Grapevine | GCA_000003745.2 |

Run `--listGenomes` to see full index paths.

---

## Example Commands

**Standard run (trim + align + QC):**
```bash
nextflow run bixBeta/wgs -r container \
  --id 7054D \
  --sheet sample-sheet.csv \
  --genome hg38 \
  --fastp \
  --bowtie2
```

**With GC bias analysis:**
```bash
nextflow run bixBeta/wgs -r container \
  --id 7054D \
  --sheet sample-sheet.csv \
  --genome mm10 \
  --fastp \
  --bowtie2 \
  --gcbias
```

**FASTQs in a subdirectory:**
```bash
nextflow run bixBeta/wgs -r container \
  --id 7054D \
  --sheet sample-sheet.csv \
  --genome hg38 \
  --fastp \
  --bowtie2 \
  --fastqs
```

**List available genomes:**
```bash
nextflow run bixBeta/wgs -r container --listGenomes
```

---

## Containers

All images are hosted on GHCR and pulled automatically by Nextflow.

| Image | Tools |
|-------|-------|
| `ghcr.io/bixbeta/wgs-align:latest` | fastp, bowtie2, samtools, picard |
| `ghcr.io/bixbeta/wgs-qualimap:latest` | qualimap |
| `ghcr.io/bixbeta/wgs-deeptools:latest` | deeptools |
| `ghcr.io/bixbeta/wgs-quarto:latest` | R 4.4.2, Quarto 1.5.57 |
| `multiqc/multiqc:v1.32` | MultiQC |

---

## Output Structure

```
./
├── trimmed_fastqs/          # fastp trimmed reads
├── trimmed_logs/            # fastp HTML + JSON reports
├── primary_BAMS/            # bowtie2 aligned BAMs
├── DEDUP_BAMS/              # deduplicated BAMs
├── STATS/
│   ├── SAMTOOLS/            # flagstat + idxstats
│   ├── PICARD/              # MarkDuplicates metrics
│   └── QUALIMAP_RES/        # bamqc output
├── GCBias_DeepTools/        # GC bias PNGs + TXT (if --gcbias)
└── Reports/
    ├── <genome>/            # MultiQC HTML report
    └── gc_bias_report.html  # Quarto GC bias report (if --gcbias)
```

---

## Tips for BioHPC

- Run inside a **`screen` or `tmux` session** so the job persists if your SSH connection drops
- Submit via a **SLURM wrapper** to run Nextflow itself on a compute node:
  ```bash
  srun --pty --mem=8G --cpus-per-task=2 bash
  nextflow run bixBeta/wgs -r container ...
  ```
- Use `-resume` to restart from the last successful step after a failure:
  ```bash
  nextflow run bixBeta/wgs -r container ... -resume
  ```
- Nextflow work files are written to `./work/` — clean up with `nextflow clean -f` after a successful run
