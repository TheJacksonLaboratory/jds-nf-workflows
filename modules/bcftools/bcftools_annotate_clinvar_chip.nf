process BCFTOOLS_ANNOTATE_CLINVAR {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '02:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    // clonal_hema only (new module; no shared jds module is touched). Adds the ClinVar INFO fields
    // (CLNSIG, CLNREVSTAT, CLNDN, CLNDISDB, GENEINFO, ALLELEID) to the variants of one VCF.
    // No publishDir: runs per SNP/INDEL branch (subworkflows/chip_annotation.nf); only the final
    // merged, fully-annotated VCF (gatk_mergevcf_chip.nf) is a published product.
    //
    // The ClinVar file must be a bgzipped, tabix-indexed VCF with chr-prefixed contig names, like
    // this workflow's VCFs. A contig-name mismatch is SILENT (exit 0, nothing annotated), so the
    // task stops if a variant's contig is not a contig of the ClinVar file. Zero ClinVar hits is
    // otherwise normal for a small panel and is not an error.
    //
    // Lookups are restricted to the variant positions (-R sites --regions-overlap 0); without it
    // bcftools streams the whole resource.

    input:
    tuple val(sampleID), file(vcf)
    path(clinvar)
    path(clinvar_index)
    val(output_suffix)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    def baseName = vcf.getName().replaceAll(/(\.vcf|\.vcf\.gz)$/, '')
    output_name = "${baseName}_${output_suffix}.vcf"
    """
    bcftools view -Oz -o in.vcf.gz ${vcf}
    bcftools index -t in.vcf.gz

    # variant positions; with no variants one dummy position keeps the header so the new INFO
    # definitions are present and the SNP and INDEL files merge cleanly
    bcftools query -f '%CHROM\\t%POS\\n' in.vcf.gz | uniq > sites.tsv
    [ -s sites.tsv ] || printf 'chr1\\t1\\n' > sites.tsv

    # fail loudly on a contig-name mismatch (e.g. '1' vs 'chr1'), which annotates nothing and exits 0
    bcftools index -s ${clinvar} | cut -f1 | sort -u > clinvar_contigs.txt
    bcftools query -f '%CHROM\\n' in.vcf.gz | sort -u > input_contigs.txt
    if [ -s input_contigs.txt ] && [ -n "\$(comm -23 input_contigs.txt clinvar_contigs.txt)" ]; then
        echo "ERROR: input contigs not found in the ClinVar file (chr-prefix mismatch?): \$(comm -23 input_contigs.txt clinvar_contigs.txt | tr '\\n' ' ')" >&2
        exit 1
    fi

    bcftools annotate -a ${clinvar} \\
        -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/GENEINFO,INFO/ALLELEID \\
        -R sites.tsv --regions-overlap 0 in.vcf.gz -O v -o ${output_name}

    N_IN=\$(bcftools view -H in.vcf.gz | wc -l)
    N_CLN=\$(grep -v '^#' ${output_name} | grep -c 'CLNSIG=' || true)
    N_HIT=\$(bcftools query -f '%CHROM\\t%POS\\n' ${output_name} | wc -l)
    echo "variants=\$N_IN clinvar_annotated=\$N_CLN"
    [ "\$N_HIT" -eq "\$N_IN" ] || { echo "ERROR: record count changed (\$N_IN -> \$N_HIT)" >&2; exit 1; }
    """
}
