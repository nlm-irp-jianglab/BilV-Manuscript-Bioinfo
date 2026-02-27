# Reproducibility Guide: bil Gene Metagenomic Survey

## Overview

This repository contains the code, gene reference database, and metadata needed to reproduce the bil gene quantification from publicly available shotgun metagenomic datasets across 14 host species (1,210 samples, 20 BioProjects).

## Directory structure

```
publication/
├── README.md                        # This file
├── table_s1_bioproject_summary.tsv  # Table S1: per-BioProject summary and sample filtering
├── table_s2_per_sample.tsv          # Table S2: per-sample SRA runs, counts, and RPKM
├── code/
│   ├── 01_process_sample.sh         # Per-sample pipeline (QC → host removal → gene alignment)
│   ├── 02_build_count_table.py      # Aggregate per-sample counts into count matrix
│   └── 03_calculate_rpkm.py         # RPKM normalization
└── data/
    ├── bil_genes.fasta              # bil gene reference (6,527 sequences)
    ├── gene_annotation.tsv          # Gene classification (bilR/bilV, basal/insertion/short)
    └── host_index.conf              # Host genome bowtie2 index paths
```

## Supplementary tables

- **Table S1** (`table_s1_bioproject_summary.tsv`): Per-BioProject summary including species, study references, original and final sample counts, sample filtering criteria, body site, geographic location, and sequencing platform.

- **Table S2** (`table_s2_per_sample.tsv`): Per-sample data including BioSample ID, BioProject, SRA run accessions (comma-separated when multiple runs per sample), non-host read count, raw bil gene counts (bilR basal/insertion/short, bilV), and RPKM values.

## Software requirements

| Tool | Version | Purpose |
|------|---------|---------|
| fastp | 0.24.0 | Read QC and adapter trimming |
| Bowtie2 | 2.5.3 | Host read removal and gene alignment |
| SAMtools | 1.21 | BAM processing and read counting |
| SRA Toolkit | 3.3.0 | SRA data download (fasterq-dump) |
| Python | 3.12 | Count aggregation, RPKM calculation, plotting |
| pandas | 2.2 | Data analysis |
| matplotlib | 3.10 | Visualization |
| seaborn | 0.13 | Visualization |

## Reproducing the analysis

### Step 1: Build gene index

```bash
bowtie2-build --threads 8 data/bil_genes.fasta data/bil_genes
```

### Step 2: Obtain host reference genomes

Host genomes used for read removal are listed in `data/host_index.conf`. Pre-built iGenomes indexes (cattle, pig, dog, rat, chicken) are available from Illumina iGenomes. Custom indexes (sheep, cat, goat, rabbit, deer, mallard, turkey) were built from NCBI RefSeq assemblies listed in the methods. The human T2T-CHM13 v2.0 assembly was used for human samples.

### Step 3: Process each sample

For each SRA run, the pipeline performs: download → quality filtering → host read removal → gene alignment → per-gene read counting. SRA run accessions for each sample are listed in the `SRA_runs` column of Table S2.

```bash
# Example: process a single sample with host read removal
bash code/01_process_sample.sh \
    -r /path/to/host_bt2_index \
    -c /path/to/cached/fastq \
    SRR_ACCESSION \
    output_dir \
    data/bil_genes
```

### Step 4: Build count matrix

Aggregates per-run idxstats output into a single count matrix, merging technical replicates by BioSample.

```bash
python code/02_build_count_table.py \
    --results-dir output_dir \
    --run-to-sample run_to_sample.tsv \
    --ncbi-runinfo run_info_all.csv \
    --bioproject-meta bioproject_meta.tsv \
    --outdir counts
```

The sample run info can be retrieved using the command:

```
esearch -db sra -query SRRXXXXXX | efetch -format runinfo
```

### Step 5: Calculate RPKM

Normalizes raw counts to RPKM (Reads Per Kilobase per Million mapped reads) and aggregates by gene group.

```bash
python code/03_calculate_rpkm.py \
    --count-matrix counts/count_matrix.tsv \
    --gene-annotation data/gene_annotation.tsv \
    --sample-metadata counts/sample_metadata.tsv \
    --outdir counts
```

## Key data files

- **bil_genes.fasta**: 6,527 bil gene sequences (4,490 bilR + 2,037 bilV) compiled from 5,000 bacterial genome assemblies (NCBI GenBank). Sequence lengths: 762–3,072 bp (mean 1,789 bp).

- **gene_annotation.tsv**: Maps each gene to its class (bilR or bilV) and group (basal, insertion, short, or bilV).

## RPKM formula

```
RPKM = (read_count × 10^9) / (reads_after_host_removal × gene_length_bp)
```

where `reads_after_host_removal` is the number of QC-passed, non-host read pairs per sample.
