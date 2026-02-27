#!/usr/bin/env python3
"""
Calculate per-gene RPKM, then summarize by gene group at two levels.

Inputs:
  --count-matrix     count_matrix.tsv (gene, gene_length, sample1, sample2, ...)
  --gene-annotation  gene_annotation.tsv (gene, gene_class, gene_group)
  --sample-metadata  sample_metadata.tsv (must have reads_after_host_removal column)
  --outdir           output directory

Outputs:
  rpkm_per_gene.tsv           — per-gene RPKM matrix
  rpkm_level1_group.tsv       — summed RPKM by gene_class (bilR, bilV) per sample
  rpkm_level2_group.tsv       — summed RPKM by gene_group (basal, insertion, short, bilV) per sample
"""

import argparse
import csv
import os
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--count-matrix', required=True)
    parser.add_argument('--gene-annotation', required=True)
    parser.add_argument('--sample-metadata', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load gene annotation
    annot = {}
    with open(args.gene_annotation) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            annot[row['gene']] = {
                'gene_class': row['gene_class'],
                'gene_group': row['gene_group'],
            }

    # Load sample normalization factors (reads_after_host_removal)
    norm = {}
    with open(args.sample_metadata) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            bs = row['BioSample']
            val = row.get('reads_after_host_removal', '0')
            norm[bs] = int(val) if val else 0

    # Read count matrix
    with open(args.count_matrix) as f:
        header = f.readline().strip().split('\t')
        samples = header[2:]  # skip gene, gene_length

        genes = []
        lengths = []
        counts = []  # list of lists

        for line in f:
            parts = line.strip().split('\t')
            genes.append(parts[0])
            lengths.append(int(parts[1]))
            counts.append([int(x) for x in parts[2:]])

    n_genes = len(genes)
    n_samples = len(samples)
    print(f"Loaded {n_genes} genes x {n_samples} samples")

    # Check normalization coverage
    missing_norm = [s for s in samples if norm.get(s, 0) == 0]
    if missing_norm:
        print(f"Warning: {len(missing_norm)} samples with 0 reads_after_host_removal (RPKM will be 0)")

    # Calculate per-gene RPKM
    # RPKM = (count * 1e9) / (reads_after_host_removal * gene_length)
    rpkm = []
    for i in range(n_genes):
        gene_len = lengths[i]
        row = []
        for j in range(n_samples):
            n = norm.get(samples[j], 0)
            if n > 0 and gene_len > 0:
                val = (counts[i][j] * 1e9) / (n * gene_len)
            else:
                val = 0.0
            row.append(val)
        rpkm.append(row)

    # Write per-gene RPKM
    rpkm_path = os.path.join(args.outdir, 'rpkm_per_gene.tsv')
    with open(rpkm_path, 'w') as f:
        f.write('gene\tgene_length\tgene_class\tgene_group\t' + '\t'.join(samples) + '\n')
        for i in range(n_genes):
            g = genes[i]
            a = annot.get(g, {'gene_class': 'bilV', 'gene_group': 'bilV'})
            vals = '\t'.join(f'{v:.4f}' for v in rpkm[i])
            f.write(f'{g}\t{lengths[i]}\t{a["gene_class"]}\t{a["gene_group"]}\t{vals}\n')
    print(f"Wrote {rpkm_path}")

    # Sum RPKM by gene_class (Level 1: bilR, bilV)
    level1_groups = sorted(set(a['gene_class'] for a in annot.values()))
    if 'bilV' not in level1_groups:
        level1_groups.append('bilV')
    level1_groups = sorted(level1_groups)

    level1_sums = {g: [0.0] * n_samples for g in level1_groups}
    for i in range(n_genes):
        a = annot.get(genes[i], {'gene_class': 'bilV', 'gene_group': 'bilV'})
        gc = a['gene_class']
        for j in range(n_samples):
            level1_sums[gc][j] += rpkm[i][j]

    l1_path = os.path.join(args.outdir, 'rpkm_level1_group.tsv')
    with open(l1_path, 'w') as f:
        f.write('gene_class\t' + '\t'.join(samples) + '\n')
        for g in level1_groups:
            vals = '\t'.join(f'{v:.4f}' for v in level1_sums[g])
            f.write(f'{g}\t{vals}\n')
    print(f"Wrote {l1_path}")

    # Sum RPKM by gene_group (Level 2: basal, insertion, short, bilV)
    level2_groups = sorted(set(a['gene_group'] for a in annot.values()))
    if 'bilV' not in level2_groups:
        level2_groups.append('bilV')
    level2_groups = sorted(level2_groups)

    level2_sums = {g: [0.0] * n_samples for g in level2_groups}
    for i in range(n_genes):
        a = annot.get(genes[i], {'gene_class': 'bilV', 'gene_group': 'bilV'})
        gg = a['gene_group']
        for j in range(n_samples):
            level2_sums[gg][j] += rpkm[i][j]

    l2_path = os.path.join(args.outdir, 'rpkm_level2_group.tsv')
    with open(l2_path, 'w') as f:
        f.write('gene_group\t' + '\t'.join(samples) + '\n')
        for g in level2_groups:
            vals = '\t'.join(f'{v:.4f}' for v in level2_sums[g])
            f.write(f'{g}\t{vals}\n')
    print(f"Wrote {l2_path}")

    # Verify: bilR at level 1 == basal + insertion + short at level 2
    bilR_l1 = level1_sums.get('bilR', [0.0] * n_samples)
    bilR_l2_sum = [level2_sums.get('basal', [0.0]*n_samples)[j]
                   + level2_sums.get('insertion', [0.0]*n_samples)[j]
                   + level2_sums.get('short', [0.0]*n_samples)[j]
                   for j in range(n_samples)]
    max_diff = max(abs(bilR_l1[j] - bilR_l2_sum[j]) for j in range(n_samples))
    print(f"\nVerification: max |bilR_level1 - (basal+insertion+short)| = {max_diff:.6e}")
    print("Level 1 = sum of Level 2 components: {'PASS' if max_diff < 0.01 else 'FAIL'}")

    # Summary stats
    print(f"\nSummary (mean RPKM across samples):")
    for g in level1_groups:
        mean_val = sum(level1_sums[g]) / n_samples
        print(f"  {g}: {mean_val:.2f}")
    for g in level2_groups:
        mean_val = sum(level2_sums[g]) / n_samples
        print(f"  {g}: {mean_val:.2f}")


if __name__ == '__main__':
    main()
