def help(){
println '''
Parameter | Default | Description

CHIP_SOMATIC — tumor-only somatic calling for a targeted CHIP gene panel
(mouse + human), using population-based germline-resource + panel-of-normals
filtering (no matched normal, no proxy-normal substitution). See
notes/germline_resource_design.md and the chip_somatic implementation plan
for the full design reasoning.

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--csv_input | /<FILE_PATH> | CSV delimited sample sheet. Required header: sampleID,fastq_1,fastq_2. Optional: sex,lane,replicate,ind. No patient/status columns — every row is an independent tumor-only sample, there is no normal-pairing concept in this workflow.

--gen_org | human | Options: human, mouse.
--genome_build | GRCh38 | Human: GRCh38. Mouse: GRCm39.
--read_type | PE | Only 'PE' is accepted for this workflow.

--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.
--detect_adapter_for_pe | false | If true, adapter auto-detection is used for paired end data.
--trim_poly_g | false | If enabled, polyG trimming is done.
--trim_poly_x | false | If enabled, polyX trimming is done.

--ref_fa | <PATH> | Reference fasta; resolved per `--gen_org`.
--ref_fa_indices | <PATH> | Pre-compiled BWA index files; resolved per `--gen_org`.
--target_gatk | <PATH> | CHIP panel target BED/interval file used for `-L` in Mutect2. Resolved per `--gen_org`. Placeholder until the panel-scoping interval list is built (tracked separately).
--target_picard | <PATH> | Picard-format target interval list, used by PICARD_COLLECTHSMETRICS.
--bait_picard | <PATH> | Picard-format bait interval list, used by PICARD_COLLECTHSMETRICS.

--germline_resource_vcf | <PATH> | Population allele-frequency VCF for Mutect2's `--germline-resource`. Human default: existing gnomAD hg38 resource (unchanged from --workflow wes). Mouse default: an empty/near-empty B6 placeholder VCF — a deliberately *safe* default (see the implementation plan, finding 1). Override directly with a file path (e.g. the UM-HET3 germline-resource) when running a different background — there is no named `--germline-background` selector; this project's specific resource files are supplied here directly, not baked into this shared repo's config.
--pon_ref | <PATH> | Panel-of-normals VCF for Mutect2's `--panel-of-normals`. Human default: existing WES PoN (1000g_pon.hg38.vcf.gz). Mouse default: a shipped, content-free, header-only placeholder (`assets/chip_somatic/empty_pon.vcf.gz`) — no real mouse CHIP-assay PoN exists yet.
--genotype_pon | true (human) / false (mouse placeholder) | Whether Mutect2 emits (still-filtered) calls at PoN-listed sites. False for mouse while the PoN is a content-free placeholder with no real sites to genotype.
--genotype_germline | true | Whether Mutect2 emits (still-filtered) calls at germline-resource-listed sites.
--contam_ref | <PATH> or null | Common-SNP resource for GATK GetPileupSummaries/CalculateContamination. Human: existing WES resource (small_exac_common_3.hg38.vcf.gz). Mouse: null — contamination estimation is skipped entirely (see finding 3: biologically inapplicable for B6, not yet built for UM-HET3).

'''
}
