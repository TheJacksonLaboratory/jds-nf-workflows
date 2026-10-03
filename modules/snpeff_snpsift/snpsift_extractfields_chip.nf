process SNPSIFT_EXTRACTFIELDS {
    tag "$sampleID"

    cpus 1
    memory  6.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d'

    // chip_somatic counterpart of snpeff_snpsift/snpsift_extractfields.nf.
    // Field list adapted from WES's 'somatic_wes' (human) case: same FORMAT
    // fields (this is Mutect2-based like WES), plus the two INFO fields that
    // drive chip_somatic's germline-vs-somatic call (POPAF, GERMQ -- see
    // notes/gnomad_af_vaf_synthetic_test_results.md) and FILTER, since
    // FilterMutectCalls' filter reasons (germline/orientation/panel_of_normals/
    // etc.) are how a reviewer will actually triage a call. Mouse: same list
    // minus the dbNSFP_* fields (dbNSFP is human-only, see chip_somatic.config).

    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern:"*.txt.gz", mode:'copy'

    input:
    tuple val(sampleID), path(vcf)
    val(gen_org)

    output:
    tuple val(sampleID), path("*.txt.gz"), emit: txt

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    // NOTE (fixed 2026-09-30, found by a real run): SnpSift's extractFields
    // expression language does NOT understand bcftools' "INFO/FIELD" prefix
    // syntax -- it parses "INFO/POPAF" as the division "INFO / POPAF" and
    // crashes with "INFO field 'INFO' not found in VCF header". Top-level
    // INFO fields are referenced by their bare name, same as WES's own field
    // lists do (see snpeff_snpsift/snpsift_extractfields.nf) -- verified by
    // rerunning after this fix.
    base_fields = 'CHROM POS ID REF ALT QUAL FILTER POPAF GERMQ "GEN[*].AF" "GEN[*].DP" "GEN[*].AD" "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE"'
    dbnsfp_fields = ' "dbNSFP_SIFT_score" "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_MutationAssessor_score" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_1000Gp3_AF" "dbNSFP_1000Gp3_AFR_AF" "dbNSFP_1000Gp3_EUR_AF" "dbNSFP_1000Gp3_AMR_AF" "dbNSFP_1000Gp3_EAS_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF"'
    fields = gen_org == 'mouse' ? base_fields : base_fields + dbnsfp_fields

    // Write to a plain file first, THEN gzip as a separate command -- not
    // piped directly (`java ... | gzip -c > out.gz`). A pipe's exit status is
    // the LAST command's by default (no `pipefail` in this process' shell),
    // so a SnpSift crash upstream was previously masked by gzip successfully
    // compressing its empty stdin into a valid-looking, empty .txt.gz (this
    // is exactly how the INFO/POPAF bug above went undetected by Nextflow:
    // the task reported SUCCESS with a real but empty output file).
    """
    java -Xmx${my_mem}G -jar /opt/snpEff/SnpSift.jar \
    extractFields ${vcf} ${fields} \
    > ${sampleID}_annotation_table.txt
    gzip ${sampleID}_annotation_table.txt
    """
}
