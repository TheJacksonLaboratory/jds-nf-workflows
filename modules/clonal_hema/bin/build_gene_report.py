#!/usr/bin/env python3
"""Build the clonal_hema cross-sample gene/variant HTML report.

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
import re
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
    ("ANN[*].HGVS_P", "Protein Change"),
    ("ANN[*].HGVS_C", "cDNA Change"),
    ("dbNSFP_SIFT_pred", "SIFT Prediction"),
    ("FILTER", "Filter"),
    ("__dbsnp__", "dbSNP"),
    ("__cosmic__", "COSMIC"),
    ("POPAF", "Pop. Frequency (gnomAD, -log10)"),
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
              "MODERATE, LOW or MODIFIER.",
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
              "value names the filter(s) that did, joined by ';'.",
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
        "The letter belongs to the transcript shown in the row; click +N to "
        "see each transcript's own letter. dbNSFP uses Ensembl transcript "
        "IDs, so each transcript is matched to dbNSFP's transcripts of the "
        "same gene at the same protein position (aa = amino-acid number, "
        "e.g. p.Pro95His is aa 95; if those disagree, both letters are "
        "shown).\n"
        "Hover over the letter to see dbNSFP's call and score for each "
        "Ensembl transcript that has one.\n"
        "Blank = the Consequence is not missense, no prediction is "
        "available, or a mouse sample (human only).",
}

# Label of the large expand button shown under the description of every metric that has
# a METRIC_DETAILS entry (the small arrow beside the metric name toggles the same detail).
METRIC_EXPAND_LABELS = {
    "Filter": "Expand for what each filter means",
    "Impact": "Expand for what each impact level means",
    "Transcript": "Expand for why one variant has several transcripts",
    "cDNA Change": "Expand for how to read c. and n. notation",
    "SIFT Prediction": "Expand for how to read the scores and the aa notation",
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
    # SIFT score/call definitions follow the dbNSFP 4.2a readme (SIFT_score, SIFT_pred).
    "SIFT Prediction": [
        "Score : from 0 to 1. The smaller the score, the more likely the "
        "change is damaging. dbNSFP calls a change D (damaging) when the "
        "score is below 0.05 and T (tolerated) when it is 0.05 or higher.",
        "Reading the scale : 0 ... D ... 0.05 ... T ... 1. Below 0.05, the "
        "closer to 0 the stronger the damaging call. From 0.05 upward, the "
        "closer to 1 the more clearly tolerated. A T just above the cutoff "
        "(for example 0.06 to 0.10) is only marginally tolerated and should "
        "be read with caution. There is no official borderline category, so "
        "scores close to 0.05 on either side are the least clear-cut.",
        "aa : the amino-acid number, i.e. the position along the protein, "
        "the same number as in the Protein Change column. Example: "
        "p.Pro95His means the proline at position 95 is replaced by a "
        "histidine, so the hover shows aa 95.",
        "Example (hover) : \"ENST... (aa 95): D (score 0.001)\" is far "
        "below the 0.05 cutoff, a strong damaging call. \"ENST... (aa 95): "
        "T (score 0.084)\" is T, but close to the cutoff, so only marginally "
        "tolerated. Both refer to amino acid 95; the two transcripts "
        "disagree.",
        "Each score is computed for one transcript's protein, so transcripts "
        "of the same gene can disagree. A mixed result such as D, T is "
        "uncertain, not a verdict.",
        "\"+N more Ensembl transcripts with no prediction\" : the hover lists "
        "only transcripts with a call; those for which dbNSFP holds no SIFT "
        "value are counted, not listed.",
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


def sift_for_line(r):
    """SIFT letter(s) and hover text for one SnpEff transcript line.

    dbNSFP lists one prediction per *Ensembl* transcript on the VCF record
    (SnpSift writes the lists comma-separated and index-aligned). The report's
    transcripts are RefSeq, and RefSeq<->Ensembl IDs only match 1:1 for MANE
    Select, so IDs are not compared: a line gets the letter(s) of the dbNSFP
    transcripts of the SAME GENE at the SAME PROTEIN POSITION. Restricting to the
    gene also stops a neighbouring gene's prediction being attached to this one."""
    def lst(key):
        v = r.get(key, "")
        return v.split(",") if v else []
    genes, aapos = lst("dbNSFP_genename"), lst("dbNSFP_aapos")
    pred, score, ids = lst("dbNSFP_SIFT_pred"), lst("dbNSFP_SIFT_score"), lst("dbNSFP_Ensembl_transcriptid")
    n = len(genes)
    if not n or not (len(aapos) == len(pred) == len(score) == len(ids) == n):
        return "", ""
    gene = r.get("ANN[*].GENE", "").upper()
    mine = [i for i in range(n) if genes[i].upper() == gene]
    m = re.match(r"p\.[A-Za-z]{3}(\d+)", r.get("ANN[*].HGVS_P", ""))
    letters = []
    if m:
        letters = list(dict.fromkeys(
            pred[i] for i in mine if aapos[i] == m.group(1) and pred[i] not in ("", ".")))
    # Only transcripts that carry a call are listed; those without one are just counted,
    # so the hover stays short but nothing is silently hidden.
    called = [i for i in mine if pred[i] not in ("", ".")]
    tip_lines = [f"{ids[i]} (aa {aapos[i]}): {pred[i]} (score {score[i]})" for i in called]
    n_none = len(mine) - len(called)
    if n_none and called:
        tip_lines.append(f"+{n_none} more Ensembl transcript{'s' if n_none != 1 else ''} with no prediction")
    tip = ("dbNSFP SIFT calls by Ensembl transcript\n"
           "(Ensembl's own list, not the RefSeq transcripts in this row)\n"
           + "\n".join(tip_lines)) if called else ""
    return ", ".join(letters), tip


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
               + ["_ntx", "_search", "_children", "_chrom_rank", "_sift_tip"])
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
        # SIFT applies to missense changes only, and is matched per transcript line (see
        # sift_for_line); the top-level cell shows the letter of the transcript displayed.
        sift_letters, sift_tip = sift_for_line(rep)
        if "missense" not in rep.get("ANN[*].EFFECT", ""):
            sift_letters, sift_tip = "", ""
        out_row[sift_idx] = sift_letters
        children = [[r.get("ANN[*].FEATUREID", ""), r.get("ANN[*].EFFECT", ""),
                     r.get("ANN[*].IMPACT", ""), r.get("ANN[*].HGVS_P", ""),
                     r.get("ANN[*].HGVS_C", ""),
                     sift_for_line(r)[0] if "missense" in r.get("ANN[*].EFFECT", "") else ""]
                    for r in ranked]
        out_row.append(len(members))
        out_row.append(" ".join(" ".join(c) for c in children))
        out_row.append(json.dumps(children))
        out_row.append(chrom_rank(out_row[chrom_idx]))
        out_row.append(sift_tip)
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
        desc = html.escape(METRIC_DESCRIPTIONS.get(disp, "")).replace("\n", "<br>")
        link = METRIC_LINKS.get(disp)
        if link:
            desc += (f' <a href="{html.escape(link[1])}" target="_blank" '
                     f'rel="noopener">{html.escape(link[0])}</a>')
        detail = METRIC_DETAILS.get(disp)
        if not detail:
            return f'<tr><td>{html.escape(disp)}</td><td>{desc}</td></tr>'
        slug = "".join(c if c.isalnum() else "-" for c in disp.lower()).strip("-")
        target = f"detail-{slug}"
        label = html.escape(METRIC_EXPAND_LABELS.get(disp, "Expand for details"))
        # Icon-only button at the right end of the description block; the words stay as its
        # accessible name (aria-label) for screen readers.
        desc = (f'<div class="metric-desc-wrap"><div class="metric-desc-text">{desc}</div>'
                f'<button type="button" class="metric-expand-btn" data-target="{target}" '
                f'aria-expanded="false" aria-label="{label}">'
                f'<svg class="chevron" viewBox="0 0 14 14" aria-hidden="true">'
                f'<path d="M1.5 4.2h11L7 10.6z"/></svg></button></div>')
        return (
            f'<tr class="metric-row-expandable">'
            f'<td>{html.escape(disp)}</td>'
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
    --expand-blue: #1e40af;
    --accent-soft: rgba(58,102,219,0.12);
    --text: #1b2030;
    --text-muted: #636b7d;
    --border: #e3e6ee;
    --switch-off: #b8bfce;
    --hover-bg: rgba(0,0,0,0.7);
    --hover-text: #ffffff;
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
      --expand-blue: #7599f0;
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
  .report-logo {{ width: 77px; height: 77px; flex: 0 0 auto; }}
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
  .dl-btn {{
    display: inline-flex; align-items: center; gap: 7px; margin-left: 16px; padding: 5px 12px;
    font: inherit; font-size: 0.85rem; color: var(--accent); background: transparent;
    border: 1px solid var(--accent); border-radius: 8px; cursor: pointer; white-space: nowrap;
  }}
  .dl-btn:hover {{ background: var(--accent-soft); }}
  .dl-btn .chevron {{ font-size: 0.8em; margin-left: -2px; }}
  .dl-wrap {{ position: relative; display: inline-block; }}
  .dl-menu {{
    position: absolute; right: 0; top: calc(100% + 6px); z-index: 20; min-width: 290px;
    background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
    box-shadow: 0 8px 24px rgba(0,0,0,0.18); padding: 6px;
  }}
  .dl-menu[hidden] {{ display: none; }}
  .dl-menu button {{
    display: flex; flex-direction: column; align-items: flex-start; gap: 2px; width: 100%;
    background: none; border: none; border-radius: 6px; padding: 8px 10px; text-align: left;
    font: inherit; color: var(--text); cursor: pointer;
  }}
  .dl-menu button:hover, .dl-menu button:focus-visible {{ background: var(--accent-soft); outline: none; }}
  .dl-menu button span {{ font-size: 0.78rem; color: var(--text-muted); }}
  .dl-btn:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
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
  .sift-tip {{ border-bottom: 1px dotted var(--text-muted); outline-offset: 2px; }}
  /* Empty cells show a dot. Pure CSS on the empty <td>: the cell's data stays empty, so sorting,
     searching and the CSV downloads are unaffected, and the dot is not selectable/copyable text. */
  #report td:empty::before, .tx-table td:empty::before {{ content: '\\00b7'; color: var(--text-muted); font-weight: 700; }}
  /* Instant hover box (Plotly-style): frosted-glass, see-through, pointer notch, no delay. */
  .hover-box {{
    display: none; position: fixed; z-index: 1000; pointer-events: none; opacity: 0;
    max-width: 380px; padding: 9px 12px; border-radius: 8px; font-size: 0.8rem; line-height: 1.45;
    background: var(--hover-bg); color: var(--hover-text); border: none;
    text-shadow: 0 0 4px rgba(165,243,252,0.65);
    box-shadow: 0 6px 24px rgba(0,0,0,0.28);
    -webkit-backdrop-filter: blur(8px) saturate(140%); backdrop-filter: blur(8px) saturate(140%);
    transition: opacity 0.05s linear;
  }}
  .hover-box.on {{ display: block; }}
  .hover-box::after {{
    content: ''; position: absolute; left: var(--arrow-x, 50%); transform: translateX(-50%);
    border: 6px solid transparent;
  }}
  .hover-box.above::after {{ top: 100%; border-top-color: var(--hover-bg); }}
  .hover-box.below::after {{ bottom: 100%; border-bottom-color: var(--hover-bg); }}
  .hover-title {{ font-weight: 600; }}
  .hover-note {{ color: rgba(207,250,254,0.8); font-size: 0.92em; margin-bottom: 4px; }}
  .hover-row {{ white-space: nowrap; font-variant-numeric: tabular-nums; }}
  @media (prefers-reduced-motion: reduce) {{ .hover-box {{ transition: none; }} }}
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
  .metric-desc-wrap {{ display: flex; align-items: flex-end; justify-content: space-between; gap: 16px; }}
  .metric-desc-text {{ flex: 1; min-width: 0; }}
  /* Box-less disclosure triangle. The button is an invisible 28px hit area (WCAG 2.5.8 asks for
     at least 24px); only the small triangle is drawn, with a soft shape-following drop-shadow
     (box-shadow would draw a rectangle). */
  .metric-expand-btn {{
    display: inline-flex; align-items: center; justify-content: center; flex: 0 0 auto;
    width: 28px; height: 28px; padding: 0; color: var(--expand-blue);
    background: none; border: none; border-radius: 50%; cursor: pointer;
  }}
  .metric-expand-btn .chevron {{
    width: 14px; height: 14px; fill: currentColor;
    filter: drop-shadow(0 2px 3px rgba(0,0,0,0.3));
    filter: drop-shadow(0 2px 3px color-mix(in srgb, currentColor 50%, transparent));
    transition: transform 0.15s ease, filter 0.15s ease;
  }}
  .metric-expand-btn:hover .chevron {{
    transform: translateY(1px) scale(1.12);
    filter: drop-shadow(0 3px 4px rgba(0,0,0,0.38));
    filter: drop-shadow(0 3px 4px color-mix(in srgb, currentColor 60%, transparent));
  }}
  .metric-expand-btn:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 0; }}
  .metric-expand-btn[aria-expanded="true"] .chevron {{ transform: rotate(180deg); }}
  .metric-expand-btn[aria-expanded="true"]:hover .chevron {{ transform: rotate(180deg) translateY(1px) scale(1.12); }}
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
    <p class="metrics-sub">What each column in the table above means. A dot (&middot;) in a cell means there is no value.</p>
    <table class="metrics">{metric_rows}</table>
  </div>

  <p class="report-footer">Generated by JDS Nextflow pipeline <code>--workflow clonal_hema</code></p>
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
const idxSift = columns.indexOf("SIFT Prediction");
const idxSiftTip = columns.indexOf("_sift_tip");

function esc(s) {{
  return String(s).replace(/[&<>"']/g, function(c) {{
    return {{'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}}[c];
  }});
}}

function transcriptTable(children) {{
  let h = '<div class="tx-detail"><table class="tx-table"><thead><tr><th>Transcript</th>' +
          '<th>Consequence</th><th>Impact</th><th>Protein Change</th><th>cDNA Change</th><th>SIFT</th></tr></thead><tbody>';
  children.forEach(function(c) {{
    h += '<tr><td>' + esc(c[0]) + '</td><td>' + esc(c[1]) + '</td><td class="impact-' + esc(c[2]) + '">' +
         esc(c[2]) + '</td><td>' + esc(c[3]) + '</td><td>' + esc(c[4]) + '</td><td>' + esc(c[5]) + '</td></tr>';
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
      {{ targets: idxSiftTip, visible: false, searchable: false, orderable: false }},
      {{ targets: idxSift, render: function(d, type, row) {{
          if (type !== 'display' || !d) return d;
          const tip = row[idxSiftTip];
          return tip ? '<span class="sift-tip" tabindex="0" data-tip="' + esc(tip) + '">' + esc(d) + '</span>' : esc(d);
        }} }},
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

  let currentSample = '';
  function selectSample(name) {{
    currentSample = name;
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

  // Download menu, right of the entries-per-page menu. Both files cover the rows currently in
  // view (selected sample + any search), every page, in the table's column order.
  //  - variants: one row per variant and gene, exactly as the table shows it.
  //  - transcripts: one row per transcript line; the variant-level columns (position, VAF, reads,
  //    filter, dbSNP, COSMIC, ...) are repeated on each, and the transcript-level columns
  //    (Transcript, Consequence, Impact, Protein/cDNA Change, SIFT) hold that transcript's own values.
  const idxConsequence = columns.indexOf('Consequence');
  const idxProtein = columns.indexOf('Protein Change');
  const idxCdna = columns.indexOf('cDNA Change');
  const csvCell = function(v) {{
    const s = String(v === null || v === undefined ? '' : v).replace(/;/g, '; ');
    return /[",\\r\\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }};
  function downloadCsv(kind) {{
    const keep = [];
    columns.forEach(function(c, i) {{ if (c !== 'Sample' && c.charAt(0) !== '_') keep.push(i); }});
    const head = keep.map(function(i) {{ return columns[i]; }});
    head.push('SIFT by Ensembl transcript (dbNSFP)');
    if (kind === 'transcripts') head.push('Shown in table row');
    const rows = [head.map(csvCell).join(',')];
    // The hover text without its two heading lines: one "ENST (aa N): call (score)" entry per transcript.
    const siftDetail = function(r) {{ return String(r[idxSiftTip] || '').split('\\n').slice(2).join(' | '); }};
    table.rows({{ search: 'applied', order: 'applied' }}).data().each(function(r) {{
      if (kind === 'variants') {{
        rows.push(keep.map(function(i) {{ return r[i]; }}).concat([siftDetail(r)]).map(csvCell).join(','));
        return;
      }}
      JSON.parse(r[idxChildren]).forEach(function(c, n) {{
        const own = {{}};
        own[idxTranscript] = c[0]; own[idxConsequence] = c[1]; own[idxImpact] = c[2];
        own[idxProtein] = c[3]; own[idxCdna] = c[4]; own[idxSift] = c[5];
        const vals = keep.map(function(i) {{ return (i in own) ? own[i] : r[i]; }});
        vals.push(String(c[1]).indexOf('missense') >= 0 ? siftDetail(r) : '');   // like the SIFT letter: missense lines only
        vals.push(n === 0 ? 'yes' : '');   // the first transcript is the one the table row shows
        rows.push(vals.map(csvCell).join(','));
      }});
    }});
    const blob = new Blob(['\\ufeff' + rows.join('\\r\\n') + '\\r\\n'], {{ type: 'text/csv;charset=utf-8' }});
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = (currentSample || 'CHIP-Panel') + '_' + kind + '.csv';
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(function() {{ URL.revokeObjectURL(a.href); }}, 1000);
  }}

  const $dl = $('<div class="dl-wrap">' +
      '<button type="button" class="dl-btn" id="download-menu-btn" aria-haspopup="menu" aria-expanded="false">' +
      '<svg viewBox="0 0 16 16" width="14" height="14" aria-hidden="true"><path d="M8 1.5v8m0 0L4.5 6.5M8 9.5l3.5-3M2 11.5v2h12v-2" ' +
      'fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/></svg>' +
      '<span>Download CSV</span><span class="chevron" aria-hidden="true">&#9662;</span></button>' +
      '<div class="dl-menu" role="menu" hidden>' +
      '<button type="button" role="menuitem" data-kind="variants"><strong>Variants</strong>' +
      '<span>One row per variant and gene, as shown in the table</span></button>' +
      '<button type="button" role="menuitem" data-kind="transcripts"><strong>All transcripts</strong>' +
      '<span>One row per transcript; variant columns repeated on each row</span></button>' +
      '</div></div>');
  if ($topRow.children().length >= 2) {{ $topRow.children().last().append($dl); }}
  else {{ $container.prepend($dl); }}
  const closeDlMenu = function() {{
    $dl.find('.dl-menu').attr('hidden', true);
    $dl.find('#download-menu-btn').attr('aria-expanded', 'false');
  }};
  $dl.find('#download-menu-btn').on('click', function(e) {{
    e.stopPropagation();
    const open = $dl.find('.dl-menu').is('[hidden]');
    if (open) {{ $dl.find('.dl-menu').removeAttr('hidden'); $(this).attr('aria-expanded', 'true'); }}
    else {{ closeDlMenu(); }}
  }});
  $dl.find('.dl-menu button').on('click', function() {{
    downloadCsv($(this).data('kind'));
    closeDlMenu();
  }});
  $(document).on('click', closeDlMenu).on('keydown', function(e) {{ if (e.key === 'Escape') closeDlMenu(); }});
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

  // Instant hover box for the SIFT cell (replaces the browser's delayed title tooltip).
  const $hover = $('<div class="hover-box" role="tooltip"></div>').appendTo('body');
  function showHover(el) {{
    const lines = String(el.getAttribute('data-tip') || '').split('\\n');
    if (!lines[0]) return;
    let h = '<div class="hover-title">' + esc(lines[0]) + '</div>';
    if (lines[1]) h += '<div class="hover-note">' + esc(lines[1]) + '</div>';
    lines.slice(2).forEach(function(l) {{
      h += '<div class="' + (l.charAt(0) === '+' ? 'hover-note' : 'hover-row') + '">' + esc(l) + '</div>';
    }});
    $hover.html(h).addClass('on').css({{ left: 0, top: 0 }});
    const r = el.getBoundingClientRect(), w = $hover.outerWidth(), ht = $hover.outerHeight();
    const left = Math.min(Math.max(8, r.left + r.width / 2 - w / 2), window.innerWidth - w - 8);
    const above = (r.top - ht - 12 >= 8);
    const top = above ? r.top - ht - 12 : r.bottom + 12;
    // the notch always points at the letter, even when the box is pushed sideways to stay on screen
    const arrowX = Math.min(Math.max(14, r.left + r.width / 2 - left), w - 14);
    $hover.toggleClass('above', above).toggleClass('below', !above)
          .css({{ left: left, top: top, opacity: 1, '--arrow-x': arrowX + 'px' }});
  }}
  function hideHover() {{ $hover.removeClass('on').css('opacity', 0); }}
  $('#report tbody').on('mouseenter focusin', '.sift-tip', function() {{ showHover(this); }})
                    .on('mouseleave focusout', '.sift-tip', hideHover);
  $('#table-card').on('scroll', hideHover);
  $(window).on('scroll', hideHover);

  // Redraws (sort/filter/paging) rebuild the count buttons collapsed; keep the child rows in step.
  table.on('draw.dt', hideHover);
  table.on('draw.dt', function() {{
    table.rows().every(function() {{ if (this.child.isShown()) this.child.hide(); }});
  }});

  // The expand button at the right end of a description opens/closes that metric's detail row.
  $('.metric-expand-btn').on('click', function() {{
    const target = this.dataset.target;
    const open = $('#' + target).hasClass('open');
    $('#' + target).toggleClass('open', !open);
    $(this).attr('aria-expanded', open ? 'false' : 'true');
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
