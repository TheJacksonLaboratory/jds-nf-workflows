process BCFTOOLS_MIN_ALT_MOLECULES {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    // clonal_hema only (new module; no shared jds module is touched). Minimum alt-molecule rule on the
    // FilterMutectCalls VCF: a record whose largest alt allelic depth (FORMAT AD, alt alleles) is below
    // params.min_alt_molecules gets FILTER min_alt_molecules. Soft filter: every record is kept; a PASS
    // record becomes min_alt_molecules, other FILTER values get it appended (-m +).
    // On fgbio consensus BAMs each read is one molecule, so AD alt = alt molecules. Mutect2 itself has no
    // minimum alt count: one consensus read with base quality 90 can reach PASS (GTBH26-MengQ-149).
    // Checked 2026-11-16 with this container: same CHROM/POS/REF/ALT + FILTER as an independent pysam
    // implementation in 48 of 48 GTBH26-MengQ-149 Duplex VCFs (min 2).
    //
    // Published as <sampleID>.mutect2_somatic.final.vcf.gz (sample prefix ends with a dot; 2026-11-16); its annotated
    // version is <sampleID>.mutect2_somatic.annotated_final.vcf.gz (gatk_mergevcf_chip.nf). The FilterMutectCalls chip
    // modules no longer publish their VCF (<sampleID>.mutect2_somatic.filtered.vcf.gz, work dir only).

    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern: "*.mutect2_somatic.final.vcf.gz*", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(tbi)

    output:
    tuple val(sampleID), file("${sampleID}.mutect2_somatic.final.vcf.gz"), file("${sampleID}.mutect2_somatic.final.vcf.gz.tbi"), emit: vcf_tbi

    script:
    def n = params.min_alt_molecules
    """
    printf '##FILTER=<ID=min_alt_molecules,Description="Alt allele supported by fewer than ${n} molecules (largest FORMAT AD alt < ${n}; --min_alt_molecules ${n})">\\n' > min_alt_header.txt
    bcftools annotate -h min_alt_header.txt ${vcf} \\
      | bcftools filter -s min_alt_molecules -m + -e 'MAX(FORMAT/AD[0:1-])<${n}' -Oz -o ${sampleID}.mutect2_somatic.final.vcf.gz
    bcftools index -t ${sampleID}.mutect2_somatic.final.vcf.gz
    """
}
