#!/usr/bin/env python3
"""Build the chip_somatic cross-sample gene/variant HTML report.

Input: every sample's <sample>_annotation_table.txt.gz (SNPSIFT_EXTRACTFIELDS
output, already verified correct -- see notes/test_summary_human.txt H15).
These files are the single source of truth for this report: no VCF is
re-parsed, no annotation is re-derived.

Layout (per user decision, 2026-10-02/03):
- One combined report across every sample in the run (not one per sample),
  navigated via a left-hand sample list rather than a "sample" table column
  -- clicking a sample filters the one shared table to that sample's rows
  and updates the headline above the table.
- Only rows with an actual variant call -- no synthetic "no call" rows for
  panel genes/intervals that had nothing.
- Column selection/renaming: raw field syntax like "GEN[*].AF" / "ANN[*].
  GENEID" is unreadable to a clinician audience. DISPLAY_COLUMNS is a curated
  subset, in display order, each mapped to a plain-language header; dropped
  fields (QUAL, GERMQ, transcript biotype/rank/cDNA-CDS-AA positions, dbNSFP's
  secondary in-silico scores/subpopulation AFs) are not deleted from the
  underlying data -- the full per-sample annotation_table.txt.gz still has
  every field. A "Column Reference" table under the main table documents what
  each shown column means (METRIC_DESCRIPTIONS).
- POPAF is shown as-is (a -log10 population-frequency score from Mutect2's
  germline resource), NOT converted to a literal allele-frequency percentage
  -- GATK's own Mutect2 module documents that the "not in resource" default AF
  it corresponds to is dynamically adjusted by the tool itself (see
  modules/gatk/gatk_mutect2_tumorOnly.nf), so a fixed 10**(-POPAF) conversion
  is not a verified transformation.

Output: one self-contained HTML file (table data embedded as JSON; sorting/
filtering via DataTables.js loaded from a CDN -- no server, no local JS
bundle, matches how a customer would open this: just double-click the file).

Usage: build_gene_report.py <out.html> <sample1>_annotation_table.txt.gz [...]
"""
import base64
import gzip
import html
import json
import sys

# (raw field name, display header), in the order columns appear in the table
# and in the Column Reference section below it. "Sample" is not a displayed
# column -- it is kept in the underlying row data (index 0) purely so the
# left-hand sample nav can filter the shared table by it.
DISPLAY_COLUMNS = [
    ("CHROM", "Chrom"),
    ("POS", "Position"),
    ("REF", "Ref"),
    ("ALT", "Alt"),
    ("GEN[*].AF", "VAF (%)"),
    ("GEN[*].AD", "Ref/Alt Reads"),
    ("ANN[*].EFFECT", "Consequence"),
    ("ANN[*].GENE", "Gene"),
    ("ANN[*].FEATUREID", "Transcript"),
    ("ANN[*].IMPACT", "Impact"),
    ("FILTER", "Filter"),
    ("ANN[*].HGVS_P", "Protein Change"),
    ("ANN[*].HGVS_C", "cDNA Change"),
    ("__dbsnp__", "dbSNP"),
    ("__cosmic__", "COSMIC"),
    ("POPAF", "Pop. Frequency (gnomAD, -log10)"),
    ("dbNSFP_SIFT_pred", "SIFT Prediction"),
    ("GEN[*].DP", "Depth"),
]

METRIC_DESCRIPTIONS = {
    "Gene": "Gene this row's annotation refers to (SnpEff). Usually the gene "
            "the variant lies in; a variant near or within another gene's "
            "transcripts is also listed against that gene, with a "
            "Consequence such as upstream_gene_variant, "
            "downstream_gene_variant, 5_prime_UTR_variant, etc., so one "
            "position can appear under more than one gene.",
    "Transcript": "Transcript the Protein/cDNA change is reported against "
                  "(RefSeq for human, hg38; Ensembl for mouse, GRCm39). "
                  "Each row is one variant in one gene and shows the "
                  "transcript with the highest-impact consequence (ties: "
                  "first in SnpEff's output order). A +N button means the "
                  "variant was annotated against N more transcripts of the "
                  "same gene; click it to list them all.",
    "Chrom": "Chromosome (genome build: hg38 human / GRCm39 mouse).",
    "Position": "1-based genomic coordinate of the variant.",
    "Ref": "Reference allele.",
    "Alt": "Alternate (variant) allele.",
    "dbSNP": "rs ID(s) from dbSNP (build 151), the public catalogue of "
             "known sequence variants. A match means the variant has been "
             "catalogued before; it does not by itself say whether the "
             "variant is inherited or acquired, or whether it is harmful. "
             "Blank = not in dbSNP (human samples only).",
    "COSMIC": "COSV ID(s) from COSMIC (v95), the catalogue of mutations "
              "reported in cancer samples. A match means the variant has "
              "been reported in cancer; it does not by itself make it "
              "pathogenic or CHIP-relevant. Blank = not in COSMIC (human "
              "samples only). The COSMIC site may ask you to sign in.",
    "Consequence": "Predicted effect of the variant on the transcript "
                   "(SnpEff), given as a Sequence Ontology term (the same "
                   "vocabulary Ensembl VEP uses).",
    "Impact": "SnpEff's severity category for the Consequence: HIGH, "
              "MODERATE, LOW or MODIFIER. Expand for what each means.",
    "Protein Change": "Predicted amino-acid change, HGVS protein notation "
                      "(e.g. p.Lys57Glu). Blank for non-coding consequences.",
    "cDNA Change": "Nucleotide change in HGVS notation, relative to the "
                   "transcript (e.g. c.169A>G).",
    "VAF (%)": "Variant allele fraction: Mutect2's estimate of the "
              "percentage of reads at this position that carry the "
              "alternate allele. It closely follows Alt reads divided by "
              "Depth, but can differ when read counts are very low. Reads "
              "are UMI consensus reads (each built from the raw reads that "
              "share one UMI), not raw sequencer reads.",
    "Depth": "Total number of UMI consensus reads at this position "
             "(reference + alternate), not raw sequencing depth.",
    "Ref/Alt Reads": "Number of UMI consensus reads supporting the reference "
                     "and the alternate allele, respectively (e.g. 1134/213 "
                     "means 1,134 reference reads and 213 alternate reads).",
    "Filter": "Mutect2's quality verdict on the call: does it look like a "
              "real variant in the sample, or a likely artifact or inherited "
              "(germline) variant? PASS = no filter flagged it. Any other "
              "value names the filter(s) that did, joined by ';'. Expand for "
              "what each one means.",
    "Pop. Frequency (gnomAD, -log10)": "Population-frequency score from "
        "Mutect2's gnomAD germline resource, on a -log10 scale. A LOWER "
        "number means the variant is MORE common in the general "
        "population, which makes an inherited (germline) origin more "
        "likely; a HIGHER number means it is rare or absent from the "
        "database. A low-frequency match does not rule out CHIP: recurrent "
        "CHIP hotspots can also appear in gnomAD at low frequency (here, "
        "DNMT3A p.Arg882His scores 3.35).",
    "SIFT Prediction": "Computational (in-silico) prediction from dbNSFP of "
        "whether a missense change is damaging to the protein (D) or "
        "tolerated (T). It is a prediction, not a clinical classification. "
        "If dbNSFP's transcripts disagree, both letters are shown. Shown "
        "only on rows whose Consequence is missense; blank otherwise, when "
        "no prediction is available, and for mouse samples (human only).",
}

# Optional external reference appended to a metric's description:
# display header -> (link text, URL).
METRIC_LINKS = {
    "dbSNP": ("NCBI dbSNP", "https://www.ncbi.nlm.nih.gov/snp/"),
    "COSMIC": ("COSMIC website", "https://cancer.sanger.ac.uk/cosmic"),
    "Consequence": ("Definitions of every term (SnpEff effect table)",
                    "https://pcingola.github.io/SnpEff/snpeff/inputoutput/"
                    "#variant-annotation-details"),
}

# Optional extra detail for a metric, shown as a click-to-expand row under its
# Column Reference entry (not everything needs one -- only added where the
# short description in METRIC_DESCRIPTIONS leaves a genuinely non-obvious
# convention unexplained). Content below is checked directly against
# hgvs-nomenclature.org (2026-10-03), not recalled from memory.
METRIC_DETAILS = {
    # Filter wording follows the ##FILTER definitions in the Mutect2
    # (FilterMutectCalls) VCF header of this pipeline's own output.
    "Filter": [
        "PASS : none of Mutect2's filters flagged the call, so it is not "
        "judged to be a technical artifact or an inherited variant. This is "
        "a technical-quality verdict only: it does not say the variant is "
        "clinically significant or CHIP-relevant, and with no matched normal "
        "sample 'somatic' is inferred, not proven.",
        "germline : the evidence indicates the variant is inherited (present "
        "in the person's germline), not acquired. Judged mainly from how "
        "common the variant is in the population (gnomAD) and its allele "
        "fraction.",
        "weak_evidence : the call does not reach Mutect2's likelihood "
        "threshold, i.e. too little support to call it confidently.",
        "strand_bias : the reads supporting the variant come from one read "
        "direction only, which is more typical of an artifact than a true "
        "variant.",
        "orientation : Mutect2's orientation-bias model detected a read-"
        "orientation pattern more consistent with an artifact than a true "
        "variant.",
        "map_qual : the reads supporting the variant align noticeably worse "
        "than the reference reads, so the variant may come from "
        "misalignment.",
        "panel_of_normals : the site is on the panel-of-normals blacklist, "
        "i.e. it recurs in normal samples and is treated as background.",
        "A call can carry several filters (e.g. orientation;weak_evidence). "
        "Other filters (base_qual, clustered_events, contamination, "
        "duplicate, low_allele_frac, slippage, and more) are defined in the "
        "Mutect2 VCF header.",
    ],
    # Impact wording follows the "Impact" table of the SnpEff documentation.
    "Impact": [
        "HIGH : assumed to be disruptive to the protein, probably causing "
        "truncation, loss of function or nonsense-mediated decay (e.g. "
        "stop_gained, frameshift_variant).",
        "MODERATE : a non-disruptive change that might alter how well the "
        "protein works (e.g. missense_variant, inframe_deletion).",
        "LOW : assumed to be mostly harmless or unlikely to change protein "
        "behaviour (e.g. synonymous_variant).",
        "MODIFIER : usually non-coding variants or variants in non-coding "
        "genes, where predictions are difficult or there is no evidence of "
        "impact (e.g. downstream_gene_variant).",
    ],
    "Transcript": [
        "A gene usually has several transcripts (isoforms) built from "
        "different start sites or exons. The same genomic change is "
        "annotated against each one, and its position and effect can "
        "differ between them.",
        "Example from this assay's DNMT3A hotspot: the one genomic change "
        "is c.2645G>A, p.Arg882His on NM_022552.4 but c.2189G>A, "
        "p.Arg730His on the shorter isoform NM_001320893.1. It is the same "
        "change, numbered differently.",
        "A change can also be coding in one isoform and non-coding in "
        "another, so consequences may differ between transcripts of the "
        "same gene.",
        "A variant near or within a neighbouring gene's transcripts is "
        "listed as its own row for that gene, with a Consequence such as "
        "upstream_gene_variant, downstream_gene_variant, "
        "5_prime_UTR_variant, etc.",
    ],
    "cDNA Change": [
        "c. : position counted from the A of the transcript's start codon "
        "(ATG); negative numbers count backward into the 5' untranslated "
        "region.",
        "n. : position counted from the start of a non-coding transcript's "
        "own sequence (used when the transcript has no coding sequence); a "
        "negative n. position marks a base upstream of that start.",
        "Example c.169A>G : the 169th coding base, counted from the ATG, "
        "changes from A to G.",
        "Example n.-67C>T : 67 bases upstream of the non-coding transcript's "
        "start, C changes to T.",
        "The bases shown are on the transcript's own strand. For a gene "
        "transcribed from the minus strand they are the complement of the "
        "Ref/Alt columns, which are always the genome's plus strand "
        "(e.g. genomic T>C reads as A>G).",
    ],
}


IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}


def chrom_rank(chrom):
    """Natural chromosome order: 1..N, then X, Y, mitochondrial, then any
    other contig (a plain text sort would give chr1, chr10, chr11, ... chr2)."""
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if c.isdigit():
        return int(c)
    return {"X": 1000, "Y": 1001, "M": 1002, "MT": 1002}.get(c.upper(), 2000)


def read_table(path):
    sample = path.split("/")[-1].replace("_annotation_table.txt.gz", "")
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            vals = line.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            row["sample"] = sample
            rows.append(row)
    return header, rows


def fmt_vaf_pct(value):
    """GEN[*].AF is a 0-1 fraction; show as a percentage. Leave anything that
    isn't a single plain number (blank, multi-value) exactly as extracted
    rather than guess at it."""
    try:
        return f"{float(value) * 100:.2f}"
    except (TypeError, ValueError):
        return value


def logo_img_tag(logo_path):
    """Inline the SVG as a data URI so the report stays one self-contained
    file; its internal ids/styles can't collide with the page's own."""
    if not logo_path:
        return ""
    with open(logo_path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("ascii")
    return (f'<img class="report-logo" alt="Genome Technologies logo" '
            f'src="data:image/svg+xml;base64,{b64}">')


def main():
    args = sys.argv[1:]
    logo_path = None
    if "--logo" in args:
        i = args.index("--logo")
        logo_path = args[i + 1]
        del args[i:i + 2]
    if len(args) < 2:
        raise SystemExit("usage: build_gene_report.py [--logo <logo.svg>] <out.html> <table.txt.gz> [...]")
    out_path = args[0]
    in_paths = args[1:]

    all_rows = []
    for p in in_paths:
        _, rows = read_table(p)
        all_rows.append(rows)

    raw_fields = [raw for raw, _ in DISPLAY_COLUMNS]
    # Three hidden trailing columns: the group's transcript count, searchable
    # text of every transcript in the group, and the per-transcript details
    # for the expander on the Transcript cell.
    columns = (["Sample"] + [disp for _, disp in DISPLAY_COLUMNS]
               + ["_ntx", "_search", "_children", "_chrom_rank"])
    vaf_idx = columns.index("VAF (%)")
    reads_idx = columns.index("Ref/Alt Reads")
    sift_idx = columns.index("SIFT Prediction")
    conseq_idx = columns.index("Consequence")
    pos_idx = columns.index("Position")
    chrom_idx = columns.index("Chrom")

    # One table row per (sample, variant, gene): the SnpEff output has one line
    # per transcript, so a variant in a gene with several isoforms repeats.
    groups = {}
    for rows in all_rows:
        for row in rows:
            key = (row["sample"], row.get("CHROM", ""), row.get("POS", ""),
                   row.get("REF", ""), row.get("ALT", ""),
                   row.get("ANN[*].GENE", ""))
            groups.setdefault(key, []).append(row)

    combined = []
    sample_counts = {}
    for members in groups.values():
        # sorted() is stable: equal-impact transcripts keep SnpEff's own order.
        ranked = sorted(members, key=lambda r: -IMPACT_RANK.get(r.get("ANN[*].IMPACT", ""), -1))
        rep = dict(ranked[0])
        # The VCF ID column is dbSNP (rs...) then COSMIC (COSV...) IDs joined
        # by ';'; anything that is not an rs ID can only have come from COSMIC.
        ids = [x for x in rep.get("ID", "").split(";") if x and x != "."]
        rep["__dbsnp__"] = ";".join(x for x in ids if x.startswith("rs"))
        rep["__cosmic__"] = ";".join(x for x in ids if not x.startswith("rs"))
        out_row = [rep["sample"]] + [rep.get(c, "") for c in raw_fields]
        out_row[vaf_idx] = fmt_vaf_pct(out_row[vaf_idx])
        # SnpSift writes AD as "ref,alt"; a slash can't be mistaken for a thousands separator.
        out_row[reads_idx] = out_row[reads_idx].replace(",", "/")
        # dbNSFP gives one SIFT letter per dbNSFP transcript (e.g. "D,D,D,D", "." = none);
        # show each distinct letter once.
        out_row[sift_idx] = ", ".join(dict.fromkeys(
            x for x in out_row[sift_idx].split(",") if x and x != "."))
        # SnpSift attaches dbNSFP values to the whole VCF record, so a neighbouring gene's
        # upstream/downstream row would otherwise show the missense gene's SIFT letter.
        if "missense" not in out_row[conseq_idx]:
            out_row[sift_idx] = ""
        children = [[r.get("ANN[*].FEATUREID", ""), r.get("ANN[*].EFFECT", ""),
                     r.get("ANN[*].IMPACT", ""), r.get("ANN[*].HGVS_P", ""),
                     r.get("ANN[*].HGVS_C", "")] for r in ranked]
        out_row.append(len(members))
        out_row.append(" ".join(" ".join(c) for c in children))
        out_row.append(json.dumps(children))
        out_row.append(chrom_rank(out_row[chrom_idx]))
        try:
            out_row[pos_idx] = int(out_row[pos_idx])
        except ValueError:
            pass
        combined.append(out_row)
        sample_counts[out_row[0]] = sample_counts.get(out_row[0], 0) + 1

    samples = sorted(sample_counts)

    data_json = json.dumps(combined)
    columns_json = json.dumps(columns)
    sample_counts_json = json.dumps(sample_counts)
    samples_json = json.dumps(samples)

    n_samples = len(in_paths)
    n_variants = len(combined)
    logo_html = logo_img_tag(logo_path)

    nav_items = "".join(
        f'<button class="sample-nav-item" data-sample="{html.escape(s)}">'
        f'<span class="sample-nav-name">{html.escape(s)}</span>'
        f'<span class="count-badge">{sample_counts[s]}</span></button>'
        for s in samples
    )

    def metric_row_html(disp):
        desc = html.escape(METRIC_DESCRIPTIONS.get(disp, ""))
        link = METRIC_LINKS.get(disp)
        if link:
            desc += (f' <a href="{html.escape(link[1])}" target="_blank" '
                     f'rel="noopener">{html.escape(link[0])}</a>')
        detail = METRIC_DETAILS.get(disp)
        if not detail:
            return f'<tr><td>{html.escape(disp)}</td><td>{desc}</td></tr>'
        slug = "".join(c if c.isalnum() else "-" for c in disp.lower()).strip("-")
        target = f"detail-{slug}"
        return (
            f'<tr class="metric-row-expandable">'
            f'<td><button type="button" class="metric-term-toggle" '
            f'data-target="{target}" aria-expanded="false">'
            f'{html.escape(disp)}<span class="chevron">&#9656;</span></button></td>'
            f'<td>{desc}</td></tr>'
            f'<tr class="metric-detail-row" id="{target}">'
            f'<td></td><td><ul class="metric-detail-list">'
            + "".join(f"<li>{html.escape(item)}</li>" for item in detail)
            + "</ul></td></tr>"
        )

    metric_rows = "".join(metric_row_html(disp) for _, disp in DISPLAY_COLUMNS)

    html_out = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CHIP Panel Variant Report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Lexend:wght@500;600;700&family=IBM+Plex+Sans:wght@400;500;600&display=swap" rel="stylesheet">
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/datatables.net-dt@2/css/dataTables.dataTables.min.css">
<script src="https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/datatables.net@2/js/dataTables.min.js"></script>
<style>
  :root {{
    --bg: #f4f5f9;
    --surface: #ffffff;
    --sidebar-bg: #161a2b;
    --sidebar-text: #aeb4c8;
    --sidebar-text-active: #ffffff;
    --sidebar-border: rgba(255,255,255,0.08);
    --accent: #3a66db;
    --accent-soft: rgba(58,102,219,0.12);
    --text: #1b2030;
    --text-muted: #636b7d;
    --border: #e3e6ee;
    --switch-off: #b8bfce;
    --pass: #1e8a5f;
    --pass-bg: #e3f5ec;
    --fail: #b4471f;
    --fail-bg: #fbe9e2;
    --impact-high: #b23a3a;
    --impact-moderate: #ad7a12;
    --impact-low: #1e7a4c;
    --impact-modifier: #636b7d;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg: #11141f;
      --surface: #1a1e2d;
      --sidebar-bg: #0b0d16;
      --border: #2b3044;
      --switch-off: #3b4259;
      --text: #e7e9f2;
      --text-muted: #9099b0;
      --accent: #7599f0;
      --accent-soft: rgba(117,153,240,0.16);
      --pass-bg: rgba(30,138,95,0.18);
      --fail-bg: rgba(180,71,31,0.2);
    }}
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; display: flex; min-height: 100vh;
    font-family: 'IBM Plex Sans', -apple-system, Helvetica, Arial, sans-serif;
    background: var(--bg); color: var(--text);
  }}
  h1, h2, h3, .sidebar-title {{ font-family: 'Lexend', sans-serif; }}

  .sidebar {{
    width: 260px; flex: 0 0 260px; background: var(--sidebar-bg);
    color: var(--sidebar-text); padding: 24px 0; position: sticky; top: 0;
    height: 100vh; overflow-y: auto;
  }}
  .sidebar-header {{ padding: 0 22px 18px; border-bottom: 1px solid var(--sidebar-border); margin-bottom: 10px; }}
  .sidebar-title {{ font-weight: 700; font-size: 1.05rem; color: #fff; margin: 0; }}
  .sidebar-subtitle {{ font-size: 0.75rem; opacity: 0.75; margin-top: 2px; }}
  .sidebar-summary {{ font-size: 0.85rem; opacity: 0.85; margin-top: 12px; }}
  .sample-nav-heading {{
    margin: 18px 0 6px; padding: 0 22px; font-size: 0.8rem; font-weight: 600;
    letter-spacing: 0.08em; text-transform: uppercase; color: #fff;
  }}
  .sample-nav {{ list-style: none; margin: 0; padding: 0; }}
  .sample-nav-item {{
    display: flex; width: 100%; align-items: center; justify-content: space-between;
    gap: 10px; background: none; border: none; border-left: 3px solid transparent;
    color: inherit; text-align: left; padding: 11px 19px; font-size: 0.95rem;
    font-family: inherit; cursor: pointer;
  }}
  .sample-nav-name {{ overflow-wrap: anywhere; }}
  .sample-nav-item:hover {{ background: rgba(255,255,255,0.05); }}
  .sample-nav-item.active {{
    background: rgba(255,255,255,0.09); border-left-color: var(--accent);
    color: var(--sidebar-text-active); font-weight: 600;
  }}
  .count-badge {{
    font-size: 0.68rem; background: rgba(255,255,255,0.12); padding: 2px 8px;
    border-radius: 999px; flex: 0 0 auto;
  }}

  .main {{ flex: 1; min-width: 0; padding: 36px 40px 60px; }}
  .page-header {{ display: flex; align-items: flex-start; justify-content: space-between; gap: 20px; }}
  .page-header h1 {{ font-size: 1.55rem; line-height: 1.25; margin: 0 0 8px; }}
  .report-logo {{ width: 64px; height: 64px; flex: 0 0 auto; }}
  .report-footer {{
    margin-top: 40px; padding-top: 16px; border-top: 1px solid var(--border);
    color: var(--text-muted); font-size: 0.8rem;
  }}
  .page-desc {{ color: var(--text-muted); max-width: 740px; line-height: 1.55; margin: 0 0 30px; }}

  .sample-panel {{ display: flex; align-items: baseline; gap: 12px; margin-bottom: 14px; flex-wrap: wrap; }}
  .sample-panel h2 {{ font-size: 1.3rem; margin: 0; }}
  .sample-label {{ color: var(--text-muted); font-weight: 500; }}
  .sample-sub {{ color: var(--text-muted); font-size: 0.85rem; }}

  .table-card {{
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 12px; padding: 16px; overflow-x: auto;
  }}
  /* Optional freezing, switched by the Regular/Frozen radios (default: regular).
     A frozen header needs the card itself to scroll vertically. */
  .table-card.freeze-header {{ overflow-y: auto; max-height: 80vh; }}
  .table-card.freeze-header table.dataTable thead th {{
    position: sticky; top: 0; z-index: 3; background: var(--surface);
  }}
  .table-card.freeze-gene table.dataTable td.col-gene,
  .table-card.freeze-gene table.dataTable th.col-gene {{
    position: sticky; left: 0; background: var(--surface); border-right: 1px solid var(--border);
  }}
  .table-card.freeze-gene table.dataTable td.col-gene {{ z-index: 2; }}
  .table-card.freeze-gene table.dataTable thead th.col-gene {{ z-index: 4; }}

  /* On/off switches in the table toolbar (between Search and entries-per-page). */
  .dt-tools {{ display: flex; align-items: center; flex-wrap: wrap; gap: 6px 22px; margin: 0 16px; }}
  .tool-switch {{
    position: relative; display: inline-flex; align-items: center; gap: 8px;
    cursor: pointer; user-select: none; white-space: nowrap;
  }}
  .tool-switch input {{ position: absolute; opacity: 0; width: 1px; height: 1px; }}
  .switch-track {{
    position: relative; flex: 0 0 auto; width: 34px; height: 20px; border-radius: 999px;
    background: var(--switch-off); transition: background 0.15s ease;
  }}
  .switch-track::after {{
    content: ''; position: absolute; top: 2px; left: 2px; width: 16px; height: 16px;
    border-radius: 50%; background: #fff; box-shadow: 0 1px 2px rgba(0,0,0,0.25);
    transition: transform 0.15s ease;
  }}
  .tool-switch input:checked + .switch-track {{ background: var(--accent); }}
  .tool-switch input:checked + .switch-track::after {{ transform: translateX(14px); }}
  .tool-switch input:focus-visible + .switch-track {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
  @media (prefers-reduced-motion: reduce) {{
    .switch-track, .switch-track::after {{ transition: none; }}
  }}
  table.dataTable {{ font-size: 0.83em; width: 100% !important; }}
  table.dataTable th, table.dataTable td {{ white-space: nowrap; padding: 6px 10px !important; }}
  table.dataTable td a {{ color: var(--accent); text-decoration: none; }}
  table.dataTable td a:hover {{ text-decoration: underline; }}
  .pill {{ display: inline-block; padding: 2px 9px; border-radius: 999px; font-size: 0.82em; font-weight: 600; }}
  .pill-pass {{ color: var(--pass); background: var(--pass-bg); }}
  .pill-fail {{ color: var(--fail); background: var(--fail-bg); }}
  .impact-HIGH {{ color: var(--impact-high); font-weight: 600; }}
  .impact-MODERATE {{ color: var(--impact-moderate); font-weight: 600; }}
  .impact-LOW {{ color: var(--impact-low); }}
  .impact-MODIFIER {{ color: var(--impact-modifier); }}

  .metrics-card {{
    margin-top: 32px; background: var(--surface); border: 1px solid var(--border);
    border-radius: 12px; padding: 22px 26px; max-width: 900px;
  }}
  .metrics-card h3 {{ font-size: 1.05rem; margin: 0 0 4px; }}
  .metrics-card .metrics-sub {{ color: var(--text-muted); font-size: 0.82rem; margin: 0 0 16px; }}
  table.metrics {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; }}
  table.metrics td {{ text-align: left; padding: 9px 10px; border-bottom: 1px solid var(--border); vertical-align: top; }}
  table.metrics tr:last-child td {{ border-bottom: none; }}
  table.metrics td:first-child {{ white-space: nowrap; font-weight: 600; width: 1%; padding-right: 20px; }}
  table.metrics td:last-child {{ color: var(--text-muted); }}
  table.metrics a {{ color: var(--accent); }}
  .metric-term-toggle {{
    display: flex; align-items: center; gap: 6px; background: none; border: none;
    font: inherit; font-weight: 600; color: inherit; cursor: pointer; padding: 0;
  }}
  .metric-term-toggle .chevron {{
    font-size: 0.7em; color: var(--text-muted); display: inline-block;
    transition: transform 0.15s ease;
  }}
  .metric-term-toggle[aria-expanded="true"] .chevron {{ transform: rotate(90deg); }}
  .metric-row-expandable td {{ border-bottom: none; }}
  .metric-detail-row {{ display: none; }}
  .metric-detail-row.open {{ display: table-row; }}
  .metric-detail-row td:last-child {{ font-size: 0.82rem; padding-top: 0; line-height: 1.5; }}
  .metric-detail-list {{ margin: 0; padding-left: 18px; }}
  .metric-detail-list li {{ margin-bottom: 4px; }}

  .tx-toggle {{
    background: var(--accent-soft); color: var(--accent); border: none; cursor: pointer;
    font: inherit; font-weight: 600; padding: 2px 9px; border-radius: 999px;
  }}
  .tx-toggle .chevron {{ display: inline-block; font-size: 0.75em; margin-right: 3px; transition: transform 0.15s ease; }}
  .tx-toggle[aria-expanded="true"] .chevron {{ transform: rotate(90deg); }}
  .tx-detail {{ padding: 4px 0 10px 24px; }}
  table.tx-table {{ border-collapse: collapse; font-size: 0.95em; }}
  table.tx-table th, table.tx-table td {{
    text-align: left; padding: 4px 16px 4px 0; border-bottom: 1px solid var(--border); white-space: nowrap;
  }}
  table.tx-table th {{ color: var(--text-muted); font-weight: 600; }}

  @media (max-width: 760px) {{
    body {{ flex-direction: column; }}
    .sidebar {{ width: 100%; flex: none; position: static; height: auto; }}
    .main {{ padding: 24px 16px 40px; }}
  }}
</style>
</head>
<body>
<nav class="sidebar">
  <div class="sidebar-header">
    <p class="sidebar-title">CHIP Panel</p>
    <p class="sidebar-subtitle">Variant Report</p>
    <p class="sidebar-summary">{n_samples} sample(s) &middot; {n_variants} variant(s)</p>
  </div>
  <p class="sample-nav-heading">Select a sample</p>
  <div class="sample-nav">{nav_items}</div>
</nav>
<main class="main">
  <div class="page-header">
    <div>
      <h1>Clonal Hematopoiesis of Indeterminate Potential Panel<br>Variant Report</h1>
      <p class="page-desc">This report lists every variant called across the panel for each
      sample in this run, including calls that Mutect2's filters flagged (see the Filter
      column), with their predicted consequence, protein/DNA change, variant allele
      fraction, and population-frequency annotation.</p>
    </div>
    {logo_html}
  </div>

  <div class="sample-panel">
    <h2><span class="sample-label">Sample:</span> <span id="sample-headline"></span></h2>
    <span class="sample-sub" id="sample-subhead"></span>
  </div>

  <div class="table-card" id="table-card">
    <table id="report" class="display" style="width:100%">
    <thead><tr>{"".join(f"<th>{html.escape(c)}</th>" for c in columns)}</tr></thead>
    </table>
  </div>

  <div class="metrics-card">
    <h3>Column Reference</h3>
    <p class="metrics-sub">What each column in the table above means.</p>
    <table class="metrics">{metric_rows}</table>
  </div>

  <p class="report-footer">Generated by JDS Nextflow pipeline <code>--workflow chip_somatic</code></p>
</main>
<script>
const columns = {columns_json};
const data = {data_json};
const sampleCounts = {sample_counts_json};
const samples = {samples_json};

const idxSample = columns.indexOf("Sample");
const idxImpact = columns.indexOf("Impact");
const idxFilter = columns.indexOf("Filter");
const idxTranscript = columns.indexOf("Transcript");
const idxNtx = columns.indexOf("_ntx");
const idxSearch = columns.indexOf("_search");
const idxChrom = columns.indexOf("Chrom");
const idxPos = columns.indexOf("Position");
const idxChromRank = columns.indexOf("_chrom_rank");
const idxDbsnp = columns.indexOf("dbSNP");
const idxCosmic = columns.indexOf("COSMIC");
const idxGene = columns.indexOf("Gene");

function idLinks(d, type, urlPrefix, linkIf) {{
  if (type !== 'display' || !d) return d;
  return d.split(';').map(function(id) {{
    return linkIf(id)
      ? '<a href="' + urlPrefix + encodeURIComponent(id) + '" target="_blank" rel="noopener">' + esc(id) + '</a>'
      : esc(id);
  }}).join(', ');
}}
const idxChildren = columns.indexOf("_children");

function esc(s) {{
  return String(s).replace(/[&<>"']/g, function(c) {{
    return {{'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}}[c];
  }});
}}

function transcriptTable(children) {{
  let h = '<div class="tx-detail"><table class="tx-table"><thead><tr><th>Transcript</th>' +
          '<th>Consequence</th><th>Impact</th><th>Protein Change</th><th>cDNA Change</th></tr></thead><tbody>';
  children.forEach(function(c) {{
    h += '<tr><td>' + esc(c[0]) + '</td><td>' + esc(c[1]) + '</td><td class="impact-' + esc(c[2]) + '">' +
         esc(c[2]) + '</td><td>' + esc(c[3]) + '</td><td>' + esc(c[4]) + '</td></tr>';
  }});
  return h + '</tbody></table></div>';
}}

function escapeRegex(s) {{
  return s.replace(/[.*+?^${{}}()|[\\]\\\\]/g, '\\\\$&');
}}

$(document).ready(function() {{
  const table = $('#report').DataTable({{
    data: data,
    columns: columns.map(c => ({{ title: c }})),
    columnDefs: [
      {{ targets: idxSample, visible: false, searchable: true }},
      {{ targets: idxSearch, visible: false, searchable: true }},
      {{ targets: idxChildren, visible: false, searchable: false, orderable: false }},
      {{ targets: idxNtx, visible: false, searchable: false }},
      {{ targets: idxGene, className: 'col-gene' }},
      {{ targets: idxChromRank, visible: false, searchable: false }},
      {{ targets: idxDbsnp, render: function(d, type) {{
          return idLinks(d, type, 'https://www.ncbi.nlm.nih.gov/snp/', function(id) {{ return /^rs[0-9]+$/.test(id); }});
        }} }},
      {{ targets: idxCosmic, render: function(d, type) {{
          return idLinks(d, type, 'https://cancer.sanger.ac.uk/cosmic/search?q=', function(id) {{ return /^COSV[0-9]+$/.test(id); }});
        }} }},
      {{ targets: [idxChrom, idxPos], orderData: [idxChromRank, idxPos] }},
      {{ targets: idxTranscript, render: function(d, type, row) {{
          if (type !== 'display') return d;
          const n = row[idxNtx];
          if (n < 2) return esc(d);
          return esc(d) + ' <button type="button" class="tx-toggle" aria-expanded="false" ' +
                 'title="Show all ' + n + ' transcripts"><span class="chevron">&#9656;</span>+' + (n - 1) + '</button>';
        }} }},
      {{ targets: idxImpact, createdCell: function(td, cellData) {{
          if (cellData) $(td).addClass('impact-' + cellData);
        }} }},
      {{ targets: idxFilter, createdCell: function(td, cellData) {{
          if (!cellData) return;
          const ok = cellData === 'PASS';
          $(td).html('<span class="pill ' + (ok ? 'pill-pass' : 'pill-fail') + '">' + cellData + '</span>');
        }} }}
    ],
    pageLength: 25,
    language: {{ infoFiltered: '' }},
    layout: {{ topStart: 'search', topEnd: 'pageLength', bottomStart: 'info', bottomEnd: 'paging' }},
    order: [[ idxChrom, 'asc' ]],
    deferRender: true
  }});

  function selectSample(name) {{
    $('.sample-nav-item').removeClass('active')
      .filter('[data-sample="' + CSS.escape(name) + '"]').addClass('active');
    $('#sample-headline').text(name);
    $('#sample-subhead').text((sampleCounts[name] || 0) + ' variant(s)');
    table.column(idxSample).search('^' + escapeRegex(name) + '$', true, false).draw();
  }}

  $('.sample-nav-item').on('click', function() {{
    selectSample($(this).data('sample').toString());
  }});

  $(table.column(idxGene).header()).addClass('col-gene');

  // Toolbar switches (default off), placed between the Search box and the entries-per-page menu.
  const toolSwitch = function(id, label) {{
    return '<label class="tool-switch"><input type="checkbox" role="switch" id="' + id + '">' +
           '<span class="switch-track" aria-hidden="true"></span><span>' + label + '</span></label>';
  }};
  const $tools = $('<div class="dt-tools">' + toolSwitch('freeze-header-switch', 'Freeze header') +
                   toolSwitch('freeze-gene-switch', 'Freeze Gene column') + '</div>');
  const $container = $(table.table().container());
  const $topRow = $container.find('.dt-layout-row').first();
  if ($topRow.children().length >= 2) {{ $topRow.children().first().after($tools); }}
  else {{ $container.prepend($tools); }}
  $('#freeze-header-switch').on('change', function() {{
    $('#table-card').toggleClass('freeze-header', this.checked);
  }});
  $('#freeze-gene-switch').on('change', function() {{
    $('#table-card').toggleClass('freeze-gene', this.checked);
  }});

  if (samples.length) selectSample(samples[0]);

  $('#report tbody').on('click', '.tx-toggle', function() {{
    const row = table.row($(this).closest('tr'));
    if (row.child.isShown()) {{
      row.child.hide();
      $(this).attr('aria-expanded', 'false');
    }} else {{
      row.child(transcriptTable(JSON.parse(row.data()[idxChildren]))).show();
      $(this).attr('aria-expanded', 'true');
    }}
  }});

  // Redraws (sort/filter/paging) rebuild the count buttons collapsed; keep the child rows in step.
  table.on('draw.dt', function() {{
    table.rows().every(function() {{ if (this.child.isShown()) this.child.hide(); }});
  }});

  $('.metric-term-toggle').on('click', function() {{
    const open = $(this).attr('aria-expanded') === 'true';
    $(this).attr('aria-expanded', open ? 'false' : 'true');
    $('#' + this.dataset.target).toggleClass('open', !open);
  }});
}});
</script>
</body>
</html>
"""
    with open(out_path, "w") as f:
        f.write(html_out)
    print(f"{n_samples} samples, {n_variants} variant rows -> {out_path}")


if __name__ == "__main__":
    main()
