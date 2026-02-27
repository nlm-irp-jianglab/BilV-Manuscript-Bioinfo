#!/usr/bin/env python3
"""
Build a comprehensive count matrix + sample metadata table.

Outputs:
  <outdir>/count_matrix.tsv        — genes (rows) × samples (columns), merged read counts
  <outdir>/sample_metadata.tsv     — one row per sample with all available metadata
  <outdir>/combined_long.tsv       — long-format table (gene, sample, count, + all metadata)

Usage:
  python3 build_count_table.py \
      --results-dir <RESULTS_DIR> \
      --run-to-sample <RUN_TO_SAMPLE.tsv> \
      --ncbi-runinfo <run_info_all.csv> \
      --bioproject-meta <selected_bioprojects_2per_species.tsv> \
      [--sample-meta <avian_gut_mgx_metadata.tsv>] \
      --outdir <OUTPUT_DIR>
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path


def load_run_to_sample(path):
    """Load run → BioSample mapping (TSV, no header)."""
    mapping = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def load_ncbi_runinfo(path):
    """Load NCBI RunInfo CSV → dict keyed by Run accession."""
    info = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            run = row.get('Run', '').strip()
            if run:
                info[run] = row
    return info


def load_bioproject_meta(path):
    """Load bioproject-level metadata TSV → dict keyed by bioproject."""
    meta = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            bp = row.get('bioproject', '').strip()
            if bp:
                meta[bp] = row
    return meta


def load_sample_meta(path):
    """Load per-sample metadata TSV (avian-style) → dict keyed by BioSample."""
    meta = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            bs = row.get('BioSample', '').strip()
            if bs:
                meta[bs] = row
    return meta


def load_idxstats(path):
    """Load idxstats TSV → list of (gene, length, mapped, unmapped)."""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rows.append((parts[0], int(parts[1]), int(parts[2]), int(parts[3])))
    return rows


def load_read_stats(path):
    """Load read_stats.tsv → dict of step → read_pairs."""
    stats = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2 and parts[0] != 'step':
                stats[parts[0]] = parts[1]
    return stats


def main():
    parser = argparse.ArgumentParser(description='Build comprehensive count table with metadata')
    parser.add_argument('--results-dir', required=True, help='Directory with per-run output folders')
    parser.add_argument('--run-to-sample', required=True, help='TSV: run_accession <tab> sample_id')
    parser.add_argument('--ncbi-runinfo', required=True, help='NCBI RunInfo CSV (run_info_all.csv)')
    parser.add_argument('--bioproject-meta', required=True, help='Bioproject metadata TSV')
    parser.add_argument('--sample-meta', default=None, help='Per-sample metadata TSV (optional, e.g. avian)')
    parser.add_argument('--exclude-samples', default=None, help='File with BioSample IDs to exclude (col 1, # comments)')
    parser.add_argument('--outdir', required=True, help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load exclusion list
    excluded_samples = set()
    if args.exclude_samples:
        with open(args.exclude_samples) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    excluded_samples.add(line.split('\t')[0])
        print(f"Loaded {len(excluded_samples)} samples to exclude")

    # Load mappings
    print("Loading run-to-sample mapping...")
    run2sample = load_run_to_sample(args.run_to_sample)
    print(f"  {len(run2sample)} runs mapped")

    print("Loading NCBI RunInfo...")
    ncbi = load_ncbi_runinfo(args.ncbi_runinfo)
    print(f"  {len(ncbi)} runs in RunInfo")

    print("Loading bioproject metadata...")
    bp_meta = load_bioproject_meta(args.bioproject_meta)
    print(f"  {len(bp_meta)} bioprojects")

    sample_meta = {}
    if args.sample_meta:
        print("Loading per-sample metadata...")
        sample_meta = load_sample_meta(args.sample_meta)
        print(f"  {len(sample_meta)} samples")

    # Discover available idxstats files
    results_dir = Path(args.results_dir)
    run_files = {}
    for f in results_dir.glob('*/*.idxstats.tsv'):
        run = f.stem.replace('.idxstats', '')
        run_files[run] = f
    print(f"Found {len(run_files)} runs with idxstats output")

    # Get gene list from first file
    first_file = next(iter(run_files.values()))
    gene_info = load_idxstats(first_file)
    genes = [(g[0], g[1]) for g in gene_info]
    print(f"Gene count: {len(genes)}")

    # Group runs by sample
    sample_runs = defaultdict(list)
    skipped = 0
    for run, fpath in run_files.items():
        sample = run2sample.get(run)
        if sample is None:
            skipped += 1
            continue
        sample_runs[sample].append((run, fpath))
    if skipped:
        print(f"  Skipped {skipped} runs not in mapping")

    # Apply exclusion list
    if excluded_samples:
        before = len(sample_runs)
        for s in excluded_samples:
            sample_runs.pop(s, None)
        print(f"  Excluded {before - len(sample_runs)} samples")

    samples = sorted(sample_runs.keys())
    print(f"Merging counts for {len(samples)} samples...")

    # Merge counts per sample
    sample_counts = {}
    sample_read_stats = {}
    for sample in samples:
        counts = [0] * len(genes)
        merged_raw = 0
        merged_qc = 0
        merged_host = 0
        merged_aligned = 0
        for run, fpath in sample_runs[sample]:
            idx = load_idxstats(fpath)
            for i, row in enumerate(idx):
                counts[i] += row[2]
            # Load read stats if available
            stats_path = fpath.parent / f"{run}.read_stats.tsv"
            if stats_path.exists():
                rs = load_read_stats(stats_path)
                merged_raw += int(rs.get('raw', 0))
                merged_qc += int(rs.get('after_qc', 0))
                merged_host += int(rs.get('after_host_removal', 0))
                merged_aligned += int(rs.get('aligned_to_genes', 0))
        sample_counts[sample] = counts
        sample_read_stats[sample] = {
            'n_runs_merged': len(sample_runs[sample]),
            'total_reads_raw': merged_raw,
            'reads_after_qc': merged_qc,
            'reads_after_host_removal': merged_host,
            'reads_aligned_to_genes': merged_aligned,
        }

    # Build sample metadata
    print("Building sample metadata...")
    meta_columns = [
        'BioSample', 'BioProject', 'host_common_name', 'host_species',
        'geo_location', 'body_site', 'gut_region', 'treatment',
        'collection_date', 'platform', 'instrument',
        'study_title', 'paper', 'doi',
        'n_runs_merged', 'total_reads_raw', 'reads_after_qc',
        'reads_after_host_removal', 'reads_aligned_to_genes',
    ]

    sample_metadata = {}
    for sample in samples:
        meta = {c: '' for c in meta_columns}
        meta['BioSample'] = sample

        # Get info from first run in group
        first_run = sample_runs[sample][0][0]
        ncbi_row = ncbi.get(first_run, {})
        meta['BioProject'] = ncbi_row.get('BioProject', '')
        meta['host_species'] = ncbi_row.get('ScientificName', '')
        meta['platform'] = ncbi_row.get('Platform', '')
        meta['instrument'] = ncbi_row.get('Model', '')

        # Bioproject-level metadata
        bp = meta['BioProject']
        bp_row = bp_meta.get(bp, {})
        meta['host_common_name'] = bp_row.get('species', '')
        meta['geo_location'] = bp_row.get('geo_location', '')
        meta['body_site'] = bp_row.get('body_site', '')
        meta['study_title'] = bp_row.get('study_title', '')

        # Per-sample metadata (avian)
        sm_row = sample_meta.get(sample, {})
        if sm_row:
            meta['host_common_name'] = sm_row.get('Host_Common_Name', meta['host_common_name'])
            meta['host_species'] = sm_row.get('Host_Species', meta['host_species'])
            meta['geo_location'] = sm_row.get('Geographic_Location', meta['geo_location'])
            meta['gut_region'] = sm_row.get('Gut_Region', '')
            meta['treatment'] = sm_row.get('Treatment', '')
            meta['collection_date'] = sm_row.get('Collection_Date', '')
            meta['paper'] = sm_row.get('Paper', '')
            meta['doi'] = sm_row.get('DOI', '')

        # Read stats
        rs = sample_read_stats[sample]
        for k, v in rs.items():
            meta[k] = str(v)

        sample_metadata[sample] = meta

    # Write sample metadata
    meta_path = os.path.join(args.outdir, 'sample_metadata.tsv')
    with open(meta_path, 'w') as f:
        f.write('\t'.join(meta_columns) + '\n')
        for sample in samples:
            f.write('\t'.join(sample_metadata[sample][c] for c in meta_columns) + '\n')
    print(f"Wrote {meta_path}")

    # Write count matrix
    matrix_path = os.path.join(args.outdir, 'count_matrix.tsv')
    with open(matrix_path, 'w') as f:
        f.write('gene\tgene_length\t' + '\t'.join(samples) + '\n')
        for i, (gname, glen) in enumerate(genes):
            row = [gname, str(glen)] + [str(sample_counts[s][i]) for s in samples]
            f.write('\t'.join(row) + '\n')
    print(f"Wrote {matrix_path}")

    # Summary
    total_with_counts = sum(1 for s in samples if sum(sample_counts[s]) > 0)
    print(f"\nSummary:")
    print(f"  Samples: {len(samples)}")
    print(f"  Samples with >0 gene counts: {total_with_counts}")
    print(f"  Genes: {len(genes)}")
    print(f"  Output: {args.outdir}/")


if __name__ == '__main__':
    main()
