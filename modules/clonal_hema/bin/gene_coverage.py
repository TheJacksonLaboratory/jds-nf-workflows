#!/usr/bin/env python3
"""clonal_hema per-gene coverage and on-target summary (2026-11-10).

Two modes, both stdlib-only:

  sample   one sample: Picard CollectHsMetrics PER_BASE_COVERAGE files of the raw (pre-UMI-dedup) and the
           consensus (dedup) BAM -> <sample>.gene_coverage.tsv (one row per gene), plus the HsMetrics files ->
           <sample>.on_target.tsv (one row: raw and dedup on-target metrics).
  summary  all samples: the per-sample files -> run-level TSVs and MultiQC custom-content files
           (gene x sample heatmap of dedup mean depth, on-target table).

Coverage values are Picard's own: bases with MAPQ >= 20 and base quality >= 20, overlapping mate bases counted
once, duplicates excluded (Picard CollectHsMetrics defaults). Genes are taken from the target names
(<gene>_region_<n>; mouse targets are <mouseGene>_<HUMANGENE>_region_<n>; promoter targets <gene>_promoter).
On-target rate = ON_TARGET_BASES / PF_BASES_ALIGNED (Picard HsMetrics field definitions).
Per-read on-target rate (optional --read-on-target, from CLONAL_HEMA_READ_ON_TARGET) = primary mapped reads
overlapping >= 1 target base / all primary mapped reads (samtools view -c -F 0xB04 [-L targets]).
Inferred sex (optional --sex-check, from CLONAL_HEMA_SEX_CHECK): MAPQ >= 30 reads on autosomes / chrX / chrY of the raw
BAM. Call from chrY/chrX reads: < SEX_Y_TO_X_FEMALE -> female (XX), > SEX_Y_TO_X_MALE -> male (XY), in between ->
unclear (e.g. a male/female DNA mixture); fewer than SEX_MIN_X_READS chrX reads -> unclear. Cut-offs set on
GTBH26-MengQ-149 (mouse, Twist and TEQUILA panels, 2026-11-11): chrY/chrX 0.00001-0.00087 in 8 libraries of 3 female
mice, 0.0083-0.086 in 10 libraries of 3 male mice; not yet checked on human data.
"""
import argparse
import collections
import csv
import glob
import os
import re


def _one_label(name):
    base = re.sub(r'_region_\d+$', '', name)
    if base.endswith('_promoter'):
        return base.replace('_', ' ')
    if '_' in base:
        mouse, human = base.split('_', 1)
        return f'{mouse} ({human})'
    return base


def gene_label(target):
    """Picard merges overlapping target intervals and joins their names with '|'
    (e.g. 'Kmt2a_KMT2A_region_25|Kmt2a_KMT2A_region_26', 'Tert_promoter|Tert_TERT_region_1').
    Same gene -> that gene; different labels (promoter + exon) -> joined with ' + '."""
    labels = []
    for name in target.split('|'):
        lab = _one_label(name)
        if lab not in labels:
            labels.append(lab)
    return ' + '.join(labels)


SEX_Y_TO_X_FEMALE, SEX_Y_TO_X_MALE, SEX_MIN_X_READS = 0.002, 0.005, 1000


def sex_fields(path):
    r = read_tsv(path)[0]
    a, x, y = int(r['sex_autosome_reads']), int(r['sex_chrX_reads']), int(r['sex_chrY_reads'])
    ytx = y / x if x else None
    if x < SEX_MIN_X_READS or ytx is None:
        call = 'unclear'
    elif ytx < SEX_Y_TO_X_FEMALE:
        call = 'female'
    elif ytx > SEX_Y_TO_X_MALE:
        call = 'male'
    else:
        call = 'unclear'
    return {'inferred_sex': call, 'chrY_to_chrX': f'{ytx:.5f}' if ytx is not None else '',
            'chrY_per_million_autosomal': f'{1e6 * y / a:.1f}' if a else '',
            'chrX_to_autosome': f'{x / a:.4f}' if a else '',
            'sex_autosome_reads': a, 'sex_chrX_reads': x, 'sex_chrY_reads': y}


def read_hsmetrics(path):
    """Picard metrics file -> dict of the first metrics row (header line starts with BAIT_SET)."""
    lines = open(path).read().splitlines()
    for i, line in enumerate(lines):
        if line.startswith('BAIT_SET'):
            return dict(zip(line.split('\t'), lines[i + 1].split('\t')))
    raise SystemExit(f'no HsMetrics table in {path}')


def per_gene(per_base, thresholds):
    """Picard PER_BASE_COVERAGE (chrom pos target coverage) -> {gene: stats}."""
    genes = collections.OrderedDict()
    region_sum = collections.defaultdict(lambda: [0, 0])  # target -> [sum, n]
    with open(per_base) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        assert header[:4] == ['chrom', 'pos', 'target', 'coverage'], header
        for line in fh:
            _c, _p, target, cov = line.rstrip('\n').split('\t')
            cov = int(cov)
            g = genes.setdefault(gene_label(target), {'bp': 0, 'sum': 0, 'zero': 0, 'ge': [0] * len(thresholds),
                                                      'regions': set()})
            g['bp'] += 1; g['sum'] += cov; g['zero'] += cov == 0; g['regions'].add(target)
            for i, t in enumerate(thresholds):
                g['ge'][i] += cov >= t
            region_sum[target][0] += cov; region_sum[target][1] += 1
    for name, g in genes.items():
        means = {r: region_sum[r][0] / region_sum[r][1] for r in g['regions']}
        low = min(means, key=means.get)
        g['mean'] = g['sum'] / g['bp']
        g['min_region'], g['min_region_mean'] = low, means[low]
    return genes


def fmt(x, nd=1):
    return '' if x is None else (f'{x:.{nd}f}' if isinstance(x, float) else str(x))


ONT_FIELDS = ['PF_BASES_ALIGNED', 'ON_TARGET_BASES', 'PCT_SELECTED_BASES', 'PCT_OFF_BAIT', 'FOLD_ENRICHMENT',
              'MEAN_TARGET_COVERAGE', 'PCT_EXC_OFF_TARGET']


def on_target_row(hs):
    pf, on = float(hs['PF_BASES_ALIGNED']), float(hs['ON_TARGET_BASES'])
    return {'on_target_rate_pct': 100 * on / pf if pf else 0.0,
            'selected_bases_pct': 100 * float(hs['PCT_SELECTED_BASES']),
            'off_bait_pct': 100 * float(hs['PCT_OFF_BAIT']),
            'fold_enrichment': float(hs['FOLD_ENRICHMENT']),
            'mean_target_coverage': float(hs['MEAN_TARGET_COVERAGE']),
            'pf_bases_aligned': int(pf), 'on_target_bases': int(on)}


def mode_sample(a):
    th = [int(x) for x in a.thresholds.split(',')]
    d, r = per_gene(a.dedup_per_base, th), per_gene(a.raw_per_base, th)
    cols = ['sample', 'gene', 'regions', 'target_bp', 'dedup_mean_depth', 'dedup_min_region_mean',
            'dedup_min_region', 'dedup_pct_0x'] + [f'dedup_pct_ge_{t}x' for t in th] + \
           ['raw_mean_depth', 'raw_min_region_mean', 'raw_pct_0x'] + [f'raw_pct_ge_{t}x' for t in th] + \
           ['raw_to_dedup_ratio']
    with open(f'{a.sample}.gene_coverage.tsv', 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for gname, g in d.items():
            rg = r.get(gname)
            row = [a.sample, gname, len(g['regions']), g['bp'], fmt(g['mean']), fmt(g['min_region_mean']),
                   g['min_region'], fmt(100 * g['zero'] / g['bp'], 2)] + \
                  [fmt(100 * x / g['bp'], 2) for x in g['ge']]
            if rg:
                row += [fmt(rg['mean']), fmt(rg['min_region_mean']), fmt(100 * rg['zero'] / rg['bp'], 2)] + \
                       [fmt(100 * x / rg['bp'], 2) for x in rg['ge']] + \
                       [fmt(rg['mean'] / g['mean'], 2) if g['mean'] else '']
            else:
                row += [''] * (3 + len(th) + 1)
            out.write('\t'.join(str(x) for x in row) + '\n')
    rows = {'raw': on_target_row(read_hsmetrics(a.raw_hsmetrics)),
            'dedup': on_target_row(read_hsmetrics(a.dedup_hsmetrics))}
    if a.read_on_target:
        rc = read_tsv(a.read_on_target)[0]
        pads = [''] + sorted({k[len('reads_on_target'):].rsplit('_', 1)[0] for k in rc
                              if k.startswith('reads_on_target_pad')}, key=lambda p: int(p[4:]))
        for x in ('raw', 'dedup'):
            n_all = int(rc[f'reads_mapped_{x}'])
            if f'reads_primary_{x}' in rc:
                rows[x]['read_pairs'] = int(rc[f'reads_primary_{x}']) // 2
            rows[x]['reads_mapped'] = n_all
            for p in pads:   # '' = unpadded targets, '_pad<N>' = targets padded by N bp per side
                n_on = int(rc[f'reads_on_target{p}_{x}'])
                rows[x].update({f'read_on_target{p}_rate_pct': 100 * n_on / n_all if n_all else 0.0,
                                f'reads_on_target{p}': n_on})
    keys = list(rows['raw'])
    sex = sex_fields(a.sex_check) if a.sex_check else {}
    with open(f'{a.sample}.on_target.tsv', 'w') as out:
        out.write('\t'.join(['sample'] + [f'{k}_{x}' for x in ('raw', 'dedup') for k in keys] + list(sex)) + '\n')
        out.write('\t'.join([a.sample] + [fmt(rows[x][k], 2) for x in ('raw', 'dedup') for k in keys]
                            + [str(v) for v in sex.values()]) + '\n')


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def mode_summary(a):
    gene_rows = [r for p in sorted(a.gene_tsv) for r in read_tsv(p)]
    ont_rows = [r for p in sorted(a.on_target_tsv) for r in read_tsv(p)]
    samples = sorted({r['sample'] for r in gene_rows})
    genes = list(collections.OrderedDict.fromkeys(r['gene'] for r in gene_rows))
    # long format, all samples
    with open('gene_coverage_all_samples.tsv', 'w') as out:
        cols = list(gene_rows[0].keys())
        out.write('\t'.join(cols) + '\n')
        for r in gene_rows:
            out.write('\t'.join(r[c] for c in cols) + '\n')
    with open('on_target_all_samples.tsv', 'w') as out:
        cols = list(ont_rows[0].keys())
        out.write('\t'.join(cols) + '\n')
        for r in ont_rows:
            out.write('\t'.join(r[c] for c in cols) + '\n')
    val = {(r['sample'], r['gene']): r for r in gene_rows}
    # matrices: genes x samples (dedup mean depth, raw mean depth)
    for col, fn in (('dedup_mean_depth', 'gene_coverage_dedup_mean_matrix.tsv'),
                    ('raw_mean_depth', 'gene_coverage_raw_mean_matrix.tsv')):
        with open(fn, 'w') as out:
            out.write('\t'.join(['gene'] + samples) + '\n')
            for g in genes:
                out.write('\t'.join([g] + [val.get((s, g), {}).get(col, '') for s in samples]) + '\n')
    # MultiQC heatmap: rows = samples, columns = genes (custom content: first row = x categories)
    with open('gene_coverage_heatmap_mqc.tsv', 'w') as out:
        out.write("# id: 'clonal_hema_gene_coverage'\n# section_name: 'Per-gene coverage (after UMI deduplication)'\n"
                  "# description: 'Mean depth per gene on the consensus (deduplicated) BAM, from Picard CollectHsMetrics "
                  "per-base coverage (MAPQ >= 20, base quality >= 20, overlapping mates counted once). Full tables with raw "
                  "(pre-deduplication) depth and % of bases at thresholds: coverage/gene_coverage_all_samples.tsv.'\n"
                  "# plot_type: 'heatmap'\n# pconfig:\n#     title: 'Mean depth per gene (dedup)'\n"
                  "#     xlab: 'Gene'\n#     ylab: 'Sample'\n")
        out.write('\t'.join(['sample'] + genes) + '\n')
        for s in samples:
            out.write('\t'.join([s] + [val.get((s, g), {}).get('dedup_mean_depth', '0') or '0' for g in genes]) + '\n')
    # MultiQC table: on-target, raw vs dedup
    pad_keys = sorted({k[:-len('_raw')] for k in ont_rows[0] if k.startswith('read_on_target_pad') and k.endswith('_raw')},
                      key=lambda k: int(k[len('read_on_target_pad'):].split('_')[0]))
    count_keys = ['reads_mapped', 'reads_on_target'] + [k.replace('read_on_target', 'reads_on_target').replace('_rate_pct', '')
                                                        for k in pad_keys]
    keep = ['on_target_rate_pct', 'read_on_target_rate_pct'] + pad_keys + ['read_pairs'] + count_keys + ['selected_bases_pct', 'off_bait_pct', 'fold_enrichment', 'mean_target_coverage']
    with open('on_target_mqc.tsv', 'w') as out:
        out.write("# id: 'clonal_hema_on_target'\n# section_name: 'On-target rate (raw vs. after UMI deduplication)'\n"
                  "# description: 'Picard CollectHsMetrics on the raw aligned BAM (before UMI deduplication) and on the "
                  "consensus BAM. On-target rate = ON_TARGET_BASES / PF_BASES_ALIGNED (per base); read on-target rate = primary "
                  "mapped reads overlapping a target / all primary mapped reads (samtools, per read; _padN = targets padded by N bp per side); selected = on or near bait "
                  "(PCT_SELECTED_BASES); off-bait = PCT_OFF_BAIT. Percentages in %.'\n# plot_type: 'table'\n")
        hdr = [f'{k}_{x}' for k in keep for x in ('raw', 'dedup') if f'{k}_{x}' in ont_rows[0]] + \
              [k for k in ('inferred_sex', 'chrY_to_chrX', 'chrX_to_autosome') if k in ont_rows[0]]
        out.write('\t'.join(['sample'] + hdr) + '\n')
        for r in ont_rows:
            out.write('\t'.join([r['sample']] + [r[h] for h in hdr]) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('sample')
    s.add_argument('--sample', required=True)
    s.add_argument('--dedup-per-base', required=True)
    s.add_argument('--raw-per-base', required=True)
    s.add_argument('--dedup-hsmetrics', required=True)
    s.add_argument('--raw-hsmetrics', required=True)
    s.add_argument('--thresholds', default='100,500,1000')
    s.add_argument('--read-on-target', default=None)
    s.add_argument('--sex-check', default=None)
    m = sub.add_parser('summary')
    m.add_argument('--gene-tsv', nargs='+', required=True)
    m.add_argument('--on-target-tsv', nargs='+', required=True)
    a = p.parse_args()
    mode_sample(a) if a.mode == 'sample' else mode_summary(a)


if __name__ == '__main__':
    main()
