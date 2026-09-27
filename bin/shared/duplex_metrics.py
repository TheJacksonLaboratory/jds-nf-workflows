#!/usr/bin/env python3
"""
TAPS duplex sequencing dedup metrics — MultiQC custom_content TSV format.

Usage:
    taps_duplex_metrics.py <sample_id> <prededup.flagstat> <postdedup.flagstat> \\
        <raw_hsmetrics.txt> <duplex_yield_metrics.txt> <family_sizes.txt>

    (The last two args may be the full list of CollectDuplexSeqMetrics output files;
    the correct files are located by suffix.)

Outputs:
    {sample_id}_duplex_summary_mqc.tsv    — table: dedup summary metrics
    {sample_id}_duplex_familysize_mqc.tsv — bargraph: family size distribution
"""

import sys
import re
import csv
from pathlib import Path


def parse_flagstat_total(path):
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            m = re.match(r"^(\d+)", line)
            if m and line.endswith("primary"):
                return int(m.group(1))
    return 0


def parse_duplex_yield(path):
    """Return the row where fraction == 1 (or last row)."""
    last = None
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            last = row
            try:
                if abs(float(row.get("fraction", 0)) - 1.0) < 1e-9:
                    return row
            except (ValueError, TypeError):
                pass
    return last


def parse_hsmetrics_coverage(path):
    """Return (MEAN_TARGET_COVERAGE, MEDIAN_TARGET_COVERAGE) from a Picard
    CollectHsMetrics text file, run with a raised COVERAGE_CAP so these
    aren't silently truncated at Picard's own WGS/WES-tuned default of 200x."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    header_idx = next(i for i, l in enumerate(lines) if l.startswith("BAIT_SET"))
    header = lines[header_idx].split("\t")
    values = lines[header_idx + 1].split("\t")
    row = dict(zip(header, values))
    return row.get("MEAN_TARGET_COVERAGE", ""), row.get("MEDIAN_TARGET_COVERAGE", "")


def parse_family_sizes(path):
    sizes = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                sizes[int(row["family_size"])] = row
            except (ValueError, KeyError):
                pass
    return sizes


def find_suffix(candidates, suffix):
    for p in candidates:
        if Path(p).name.endswith(suffix):
            return p
    raise FileNotFoundError(f"No file ending with '{suffix}' in: {candidates}")


def write_tsv(path, comment_header, column_headers, data_row):
    with open(path, "w", newline="") as fh:
        for line in comment_header:
            fh.write(f"# {line}\n")
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(column_headers)
        writer.writerow(data_row)


def main():
    if len(sys.argv) < 7:
        sys.exit(
            f"Usage: {sys.argv[0]} <sample_id> <prededup.flagstat> <postdedup.flagstat> "
            f"<raw_hsmetrics.txt> <duplex_yield_metrics.txt> <family_sizes.txt>"
        )

    sample_id       = sys.argv[1]
    prededup_path   = sys.argv[2]
    postdedup_path  = sys.argv[3]
    raw_hsmetrics_path = sys.argv[4]
    metrics_files   = sys.argv[5:]

    duplex_yield_path = find_suffix(metrics_files, ".duplex_yield_metrics.txt")
    family_size_candidates = [
        f for f in metrics_files
        if Path(f).name.endswith(".family_sizes.txt")
        and not Path(f).name.endswith(".duplex_family_sizes.txt")
    ]
    if not family_size_candidates:
        raise FileNotFoundError(f"No '*.family_sizes.txt' found in: {metrics_files}")
    family_sizes_path = family_size_candidates[0]

    total_pre  = parse_flagstat_total(prededup_path)  // 2
    total_post = parse_flagstat_total(postdedup_path) // 2
    yield_row  = parse_duplex_yield(duplex_yield_path)
    size_rows  = parse_family_sizes(family_sizes_path)
    raw_mean_depth, raw_median_depth = parse_hsmetrics_coverage(raw_hsmetrics_path)

    fold_reduction = round(total_pre / total_post, 3) if total_post else ""

    ds_frac         = yield_row.get("ds_fraction_duplexes",       "") if yield_row else ""
    ds_frac_ideal   = yield_row.get("ds_fraction_duplexes_ideal", "") if yield_row else ""
    ds_families     = yield_row.get("ds_families",  "")               if yield_row else ""
    ds_duplexes     = yield_row.get("ds_duplexes",  "")               if yield_row else ""
    cs_families     = yield_row.get("cs_families",  "")               if yield_row else ""
    ss_families     = yield_row.get("ss_families",  "")               if yield_row else ""

    try:
        ds_frac_pct = round(float(ds_frac) * 100, 4) if ds_frac else ""
    except (ValueError, TypeError):
        ds_frac_pct = ""

    try:
        ds_frac_ideal_pct = round(float(ds_frac_ideal) * 100, 4) if ds_frac_ideal else ""
    except (ValueError, TypeError):
        ds_frac_ideal_pct = ""

    # ── Section 1: Duplex summary table ──────────────────────────────────────
    summary_header = [
        "id: 'duplex_summary'",
        "section_name: 'Duplex UMI Metrics'",
        "description: >",
        "  UMI-based deduplication statistics from fgbio CollectDuplexSeqMetrics.",
        "  The dedup fold-reduction reflects the compression achieved by grouping",
        "  reads sharing the same UMI into consensus families.",
        "help_text: >",
        "  **Dedup fold-reduction:** raw reads / consensus reads. Higher values indicate",
        "  more efficient UMI grouping (more duplicates collapsed). Typical range: 1.1–3×",
        "  for whole-genome duplex sequencing at 30–50x depth.",
        "",
        "  **Duplex yield %:** fraction of unique molecules where both strands were sequenced",
        "  (ds_duplexes / ds_families). Higher duplex yield gives better error correction.",
        "",
        "  **Ideal duplex yield %:** expected duplex yield under an idealized equal-sampling",
        "  model given the observed family size distribution. Compare to observed yield to",
        "  assess strand-sampling balance.",
        "plot_type: 'table'",
        "pconfig:",
        "  id: 'duplex_summary_table'",
        "  title: 'Duplex UMI Metrics'",
        "col_config:",
        "  Raw reads:",
        "    title: 'Raw reads'",
        "    description: 'Total reads before UMI deduplication (prededup flagstat)'",
        "    format: '{:,d}'",
        "  Consensus reads:",
        "    title: 'Consensus'",
        "    description: 'Reads after UMI deduplication (postdedup flagstat)'",
        "    format: '{:,d}'",
        "  Fold reduction:",
        "    title: 'Fold-red.'",
        "    description: 'Raw reads / consensus reads — dedup compression factor'",
        "    min: 1",
        "    max: 5",
        "    scale: 'Blues'",
        "  Duplex yield %:",
        "    title: 'Duplex %'",
        "    description: 'Observed duplex yield: ds_duplexes / ds_families'",
        "    suffix: '%'",
        "    min: 0",
        "    max: 5",
        "    scale: 'Greens'",
        "  Ideal duplex yield %:",
        "    title: 'Ideal duplex %'",
        "    description: 'Expected duplex yield under ideal equal-strand-sampling model'",
        "    suffix: '%'",
        "    min: 0",
        "    max: 5",
        "    scale: 'Greens'",
        "  DS families:",
        "    title: 'DS families'",
        "    description: 'Dual-strand families: source molecules with both strands captured (ds_families)'",
        "    format: '{:,d}'",
        "  True duplexes:",
        "    title: 'True dup.'",
        "    description: 'Molecules with both strands sequenced (ds_duplexes)'",
        "    format: '{:,d}'",
        "  SS families:",
        "    title: 'SS families'",
        "    description: 'Single-strand UMI family observations (ss_families)'",
        "    format: '{:,d}'",
        "  Raw depth (mean):",
        "    title: 'Raw depth (mean)'",
        "    description: 'Mean on-target depth before UMI consensus collapsing (raw sequencing depth as designed), uncapped'",
        "    suffix: 'x'",
        "  Raw depth (median):",
        "    title: 'Raw depth (median)'",
        "    description: 'Median on-target depth before UMI consensus collapsing, uncapped'",
        "    suffix: 'x'",
    ]
    write_tsv(
        f"{sample_id}_duplex_summary_mqc.tsv",
        summary_header,
        ["Sample", "Raw reads", "Consensus reads", "Fold reduction", "Duplex yield %", "Ideal duplex yield %", "DS families", "True duplexes", "SS families", "Raw depth (mean)", "Raw depth (median)"],
        [sample_id, total_pre, total_post, fold_reduction, ds_frac_pct, ds_frac_ideal_pct, ds_families, ds_duplexes, ss_families, raw_mean_depth, raw_median_depth],
    )

    # ── Section 2: Family size distribution bargraph ──────────────────────────
    def size_pct(sz):
        row = size_rows.get(sz, {})
        frac = row.get("ds_fraction", "")
        try:
            return round(float(frac) * 100, 4) if frac else 0
        except (ValueError, TypeError):
            return 0

    sz4plus_count = sum(
        int(r.get("ds_count", 0) or 0) for s, r in size_rows.items() if s >= 4
    )
    total_ds = sum(int(r.get("ds_count", 0) or 0) for r in size_rows.values())
    sz4plus_pct = round(sz4plus_count / total_ds * 100, 4) if total_ds else 0

    familysize_header = [
        "id: 'duplex_familysize'",
        "section_name: 'Duplex Family Size Distribution'",
        "description: >",
        "  Fraction of duplex (DS) consensus families at each UMI family size.",
        "  Family size reflects the number of read pairs sharing the same UMI.",
        "  Larger families indicate higher duplication rates (more PCR or sequencing depth).",
        "help_text: >",
        "  At typical whole-genome TAPS depths (30–50x), most families are size 1–2.",
        "  Very high proportions at size 1 indicate low library complexity.",
        "  Size 4+ families are characteristic of deep-coverage or highly duplicated libraries.",
        "plot_type: 'bargraph'",
        "pconfig:",
        "  id: 'duplex_familysize_plot'",
        "  title: 'DS Family Size Distribution'",
        "  ylab: '% of DS families'",
        "  cpswitch: false",
        "  stacking: 'normal'",
    ]
    write_tsv(
        f"{sample_id}_duplex_familysize_mqc.tsv",
        familysize_header,
        ["Sample", "Size 1", "Size 2", "Size 3", "Size 4+"],
        [sample_id, size_pct(1), size_pct(2), size_pct(3), sz4plus_pct],
    )


if __name__ == "__main__":
    main()
