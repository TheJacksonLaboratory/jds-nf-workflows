# clonal_hema variant report: methods and accountability

Applies to `<project_name>.clonalHema_variant_report.html` (`clonalHema_variant_report.html` without `--project_name`), built by `modules/clonal_hema/bin/build_gene_report.py`
(process `CHIP_GENE_REPORT`). Written 2026-10-03. Any change to a rule below must update this file
and repeat the verification it cites.

Every statement is labelled **Verified** (checked against a named source, file, or a command that was
run), **Reasoning** (a design judgment built on verified facts, not itself independently confirmed), or
**Open** (unresolved, not assumed either way).

## 1. What the report is and is not

- **Verified.** It lists every Mutect2 call after FilterMutectCalls, including calls carrying a filter
  (Filter column) and calls flagged `germline`. In the two check samples, 211 of 253 transcript lines
  were PASS and 38 were `germline`.
- **Verified.** No matched normal is used (tumor-only; `config/clonal_hema.config`).
- **Reasoning.** "Somatic" is therefore an inference, never proven. PASS means "no Mutect2 filter flagged
  the call". In the synthetic diagnostic sample, PASS calls exist at VAF 45-95%, which resemble
  heterozygous or homozygous germline variants. The report states this in its Column Reference.
- **Reasoning.** The report applies no clinical definition of clonal hematopoiesis (no VAF >= 2% threshold,
  no gene list). The clinical definition requires VAF >= 2% (`11899_2022_Article_662.pdf`); this assay detects lower fractions.
  Classification is left to the reviewer.

## 2. Data flow and tool versions

| Step | Tool / resource | Source of the version |
|---|---|---|
| UMI consensus | fgbio (container `community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa`) | module file; fgbio 4.0.0 observed inside the image |
| Calling, filtering | GATK 4.4.0.0 Mutect2 + FilterMutectCalls | container line in `gatk_mutect2_tumorOnly.nf` |
| Germline resource | `af-only-gnomad.hg38.vcf.gz` | `config/clonal_hema.config` |
| Annotation | SnpEff / SnpSift 5.1d, database `hg38` (human), `GRCm39.105` (mouse) | container line; `params.gen_ver` |
| dbSNP, COSMIC | dbSNP build 151, COSMIC v95 (human only) | reference file names in the config |
| dbNSFP | 4.2a, bare-chromosome copy `dbNSFP4.2a.gatk_formatted.nochr.txt.gz` | config; dbNSFP readme |
| ClinVar | release 2026-10-04, `clinvar_20261004.GRCh38.chr.vcf.gz` (human only) | `params.clinvar` in the config |
| Report | Python 3 in the pipeline's MultiQC container; DataTables (CDN tag `@2`, a floating tag: 2.3.8 on 2026-10-03) + jQuery 3.7.1 (pinned) from a CDN | module file; the loaded DataTables file header |

- **Verified.** Variant calling runs on the fgbio UMI consensus BAM (`subworkflows/hs_clonal_hema.nf`
  passes `CHIP_UMI_CONSENSUS.out.bam_bai` to calling). So "reads" and "Depth" in the report are UMI
  consensus reads, not raw sequencer reads.
- **Verified.** Depth equals ref reads + alt reads in all 253 check-sample records. The VAF column is
  Mutect2's allele-fraction estimate, which tracks alt/(ref+alt) closely (median gap 0.0012) but can
  differ at very low counts (largest gap 0.083, at 5 ref / 1 alt).
- **Verified.** SnpEff writes one line per transcript. After `vcfEffOnePerLine`, the check samples had
  253 lines for 54 variants (28 + 26), which the report turns into 69 variant-gene rows (38 + 31).
- **Verified.** SnpSift 5.1d `dbnsfp` removes `chr` from VCF contig names before the database lookup,
  and has no option to prevent it. With the `chr`-prefixed database it annotated 0 of 26 entries;
  with a bare-name copy of the same records it annotated all of them. The pipeline therefore uses the
  bare-name copy; VCFs keep their `chr` names. The original `omics_share` files are unmodified.

## 3. Rows and the transcript shown

- **Reasoning.** One row per (sample, variant, gene). The row shows the transcript with the highest
  SnpEff impact (HIGH > MODERATE > LOW > MODIFIER); ties keep SnpEff's own output order. All transcript
  lines of the group are kept in the "+N" expansion, so nothing is discarded.
- **Open.** Whether SnpEff's output order favours the canonical (e.g. MANE Select) transcript was not
  checked. For ties, the transcript shown is therefore not guaranteed to be the canonical one (example:
  MFSD11 at chr17:76736877 shows `5_prime_UTR_variant` among 23 transcripts, all MODIFIER).
- **Verified.** Impact definitions are quoted from the SnpEff documentation table. Consequence terms are
  Sequence Ontology terms (SnpEff documentation; Ensembl VEP documentation).
- **Verified.** The same variant is listed under each gene whose transcripts it falls in or near
  (example: chr17:76736877 under SRSF2, MFSD11, MIR636, METTL23; distances 0 to 2,996 bp).

## 4. SIFT prediction (the most derived column)

### Verified facts
- dbNSFP 4.2a readme: `SIFT_score` ranges 0 to 1, smaller = more likely damaging; `SIFT_pred` is D when
  the score is below 0.05, otherwise T. Values are `;`-separated, one per Ensembl protein.
- dbNSFP 4.2a readme: `Ensembl_proteinid` corresponds to `Ensembl_transcriptid`; `aapos` and the SIFT
  fields correspond to `Ensembl_proteinid`. In 15 records at 5 loci the transcript, protein, SIFT and
  aapos lists had equal length.
- The dbNSFP header has no RefSeq identifier column (only `Ensembl_transcriptid`, `Ensembl_proteinid`).
  The report's transcripts are RefSeq (`NM_`, `NR_`) for human.
- NCBI MANE: an exact RefSeq-Ensembl match (exon coordinates, CDS, sequence) is guaranteed only for one
  transcript per gene (MANE Select) plus a few MANE Plus Clinical (v1.5 human file: 19,363 and 74).
- NCBI `gene2ensembl` README: matches are "best match", not exact (>= 80% protein overlap, splice-site
  conditions); non-best matches are not reported.
- Tested on the check samples: mapping `NM_` to `ENST` (gene2ensembl accession-only or MANE) gave a SIFT
  letter for 37 of 66 unique missense transcript lines; matching by protein position gave 51 of 66
  (that comparison was made before the same-gene restriction below was added).

### The rule used (Reasoning)
For each SnpEff transcript line with a missense consequence, the report takes the dbNSFP entries of the
**same gene** (case-insensitive symbol) whose `aapos` **equals the protein position** in the line's HGVS
protein change, and shows their distinct SIFT letters (e.g. `D`, or `D, T`). Transcript IDs are never
compared. The top-level row shows the letter for the transcript displayed in that row. SIFT is blank on
rows whose consequence is not missense, because SnpSift attaches dbNSFP values to the whole VCF record
and a neighbouring gene's row would otherwise inherit another gene's value. The hover lists the dbNSFP
transcripts of that gene that carry a call (amino-acid position, call, score) and only counts those
without one ("+N more Ensembl transcripts with no prediction"), so nothing is silently dropped.

### Validation (check run, 2 synthetic samples)
- All 36 missense rows carry a letter (D 25, T 10, mixed 1). 90 of 120 missense transcript lines in the
  expansions (counted per sample) carry their own letter.
- Three variants were re-derived by hand from the dbNSFP records (GATA2 chr3:128481168, SRSF2
  chr17:76736877, GNAS chr20:58909366); every displayed letter matched.
- Against the exact MANE ID match: 22 of 23 MANE transcripts gave an identical letter. The 23rd (SRSF2
  NM_001195427.2) has no dbNSFP value of its own and shows `D, T` borrowed from other isoforms at the
  same position.

### Open
- Position matching is an approximation of isoform identity. Where several dbNSFP transcripts share the
  position (SRSF2: all five at aa 95), every transcript at that position shows the same merged letters.
- Not measured: how often dbNSFP's gene symbol differs from SnpEff's (such rows show no SIFT), and how
  often overlapping genes at one position (e.g. paralogs) are mixed.
- dbNSFP 4.2a is older than current Ensembl; some current transcripts are absent from it.
- SIFT is a computational prediction. It is not a clinical classification, and the report says so.

## 5. Display transformations

| Column | Transformation | Status |
|---|---|---|
| VAF (%) | Mutect2 `AF` x 100, 2 decimals | Verified against AD (see section 2) |
| Ref/Alt Reads | `ref,alt` shown as `ref/alt` | Presentation only (avoids a thousands-separator reading) |
| Pop. Frequency | Mutect2 `POPAF` (-log10 scale) shown as is | Not converted: the module documents that GATK adjusts the "not in resource" default AF, so a fixed 10^-POPAF is not a verified conversion |
| dbSNP / COSMIC | VCF `ID` split: `rs...` to dbSNP, everything else to COSMIC; each linked to the source site | Reasoning: only these two resources write IDs; all 69 rows split cleanly |
| ClinVar | the VCF INFO field `CLNSIG` (ClinVar's aggregate germline classification for the single variant), underscores shown as spaces; blank when absent | Verified: field definition in the VCF header and in ClinVar's README_VCF.txt; 21 of 21 test variants equal a direct ClinVar lookup with exact allele match (2026-10-16). Germline classification only: ClinVar's somatic clinical-impact fields are not used |
| Order | natural chromosome order, then position | Presentation only |
| cDNA / Protein | SnpEff `HGVS.c` / `HGVS.p` unmodified | HGVS rules checked at hgvs-nomenclature.org |

Dropped from display (kept in the per-sample `*_annotation_table.txt.gz`): QUAL, GERMQ, transcript
bookkeeping fields, other dbNSFP scores.

## 6. Limits of the evidence

- **Open.** All validation used two synthetic samples (`JUST_DIAGNOSTIC_V6`, `TITR3_M1600_F5`). The
  numbers above describe those samples; they are not performance claims for real specimens.
- **Open.** The page's JavaScript (expansion, hover, freeze switches, sorting) has no automated test and
  could not be run on the cluster; it was checked structurally and reviewed in a browser by the user.
- **Open.** `somatic_wes` points at the same chr-prefixed dbNSFP file and very likely has the same silent
  no-op; it was not run.
- **Open.** SnpEff's RefSeq versions are older than current NCBI versions, so a displayed `NM_` may have
  a newer version upstream; position matching does not depend on the version.

## 7. Reproducing and auditing a report

1. Run `nextflow run main.nf --workflow clonal_hema ...` (see `run_clonal_hema_*.slurm` in the project
   notes). The report is `<pubdir>/gene_report/<project_name>.clonalHema_variant_report.html`.
2. Every displayed value traces to `<pubdir>/<sample>/<sample>_annotation_table.txt.gz`, then to
   `<sample>_annotated_final.vcf.gz`, then to the Mutect2 filtered VCF.
3. Re-derive a SIFT cell: `tabix <dbNSFP nochr file> <chr>:<pos>-<pos>`, keep the entry with the same
   ref/alt, split `genename`, `aapos`, `SIFT_pred`, `SIFT_score`, `Ensembl_transcriptid` on `;`, keep the
   entries whose gene equals the row's gene and whose `aapos` equals the protein position.
4. Source-of-truth for what ran: the git commit of the pipeline, the run's `.nextflow.log`, and the
   process container lines listed in section 2.

## 8. Mouse runs (`--gen_org mouse`)

A run is one species, so a report never mixes human and mouse. The subworkflow passes the species,
and whether a germline resource was used to the report explicitly (`--species`,
`--germline-resource`); nothing is guessed from empty data.

### Verified
- The sample headline reads "Human (GRCh38) Sample:" or "Mouse (GRCm39) Sample:". The assay supports only
  these two builds, so the species (`--gen_org`) alone names the build; no build setting is passed to the report.
- Mouse transcripts are Ensembl (`ENSMUST...`, SnpEff database `GRCm39.105`); all transcripts in both
  validation runs below were `ENSMUST`.
- dbSNP, COSMIC, ClinVar and dbNSFP are human-only resources and are set to null for mouse in
  `config/clonal_hema.config`. In the 2026-09-29 mouse table (97 transcript lines) the ID field was empty
  on 97 of 97 lines and `POPAF` was 7.3 on 97 of 97 lines (the default for a variant absent from an empty
  germline resource).
- Recorded 2026-09-28 (UM-HET3 test): with the empty B6 placeholder every call gets `POPAF` 7.3, so
  inherited 100%-VAF differences from GRCm39 PASS as false somatic calls.

### Rules used (Reasoning)
- For mouse the report omits the columns that cannot carry data: dbSNP, COSMIC, ClinVar, SIFT Prediction and
  Pop. Frequency. The set is declared per species, not detected, so a missing human annotation can never
  be hidden by it. The Column Reference, transcript expansion and CSV downloads follow the same column set.
- If the run used the empty placeholder germline resource (detected by the file name
  `empty_germline_resource...`), the report shows a banner, the Filter text says inherited strain variants
  are not flagged, and `germline` is described as unavailable. If a real germline resource was supplied
  (for example UM-HET3), the banner is absent and `germline` is described as judged from that resource.
- The `panel_of_normals` bullet is dropped for mouse: the mouse panel of normals is a content-free
  placeholder (config comment, 2026-09-28).

### Validation (real runs, 2026-10-05, one TITR3 sample each, 28 tasks, 0 failed)
- B6 placeholder run: 22 variant-gene rows, all PASS, banner shown, 14 visible columns.
- UM-HET3 resource run: 58 rows (PASS 28, `germline` 29, `strand_bias` 1), no banner, 14 visible columns.
  The two runs used different input FASTQs (UM-HET3 test vs B6 titration), so these counts are not a
  like-for-like comparison of the two resources.
- The human report was regenerated unchanged apart from the label and the single-build Chrom text.

### Open
- No mouse equivalent of SIFT, dbSNP or COSMIC is shown. The one dbSNP-like resource found (MGP) is on
  GRCm38 and would need a liftover to GRCm39.
- No review flag for inherited-looking calls is built (decided 2026-10-05): the table already sorts and
  searches by VAF, so a flag would repeat that and needed a threshold with no data behind it. Instead the
  no-germline-resource banner tells the reader that inherited variants in an inbred strain typically read
  near 100% VAF and to sort the VAF column. In the UM-HET3 test the 6 variants at 99.9-100% were exactly the
  ones a real germline resource filtered; 12 of the 18 `germline` variants sat at about 49-60%, so a VAF
  cut-off alone cannot find heterozygous inherited variants. The real remedy is a strain germline resource
  on the calling side.
- Whether the UM-HET3 `germline` calls are correct was not re-validated here.
- Banner detection depends on the placeholder's file name; another empty resource under a different name
  would not trigger it.
- Mouse gene symbols are shown as SnpEff reports them; no human ortholog label is added (decision 2026-10-05).
