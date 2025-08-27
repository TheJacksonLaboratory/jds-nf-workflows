process GATK_FILTERMUECTCALLS {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID + '/mt_callers'}", pattern: "*exclusionFiltered.final.vcf*", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/mt_callers/stats'}", pattern: "*.filteringStats.tsv", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(stats), path(contamination)
    val(type)

    output:
    tuple val(sampleID), file("*.exclusionFiltered*.vcf.gz"), file("*.exclusionFiltered*.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(sampleID), file("*.filteringStats.tsv"), emit: stats

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    def extra_args = ""
    if (contamination.toString().toLowerCase().contains('mouse')) {
        extra_args = "--max-alt-allele-count ${params.max_allele_count}"
    } else if (contamination.toString() == 'primary') {
        extra_args = "--min-allele-fraction 0 --max-alt-allele-count 4"
    } else {
        extra_args = "--contamination-estimate `cat ${contamination}` --max-alt-allele-count ${params.max_allele_count}"
    }

    """
    mkdir -p tmp

    if [[ "${vcf}" == *.vcf.gz ]]; then
        gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" IndexFeatureFile \
        -I ${vcf}
    fi

    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" FilterMutectCalls \
        -R ${params.mt_fasta} \
        -V ${vcf} \
        --stats ${stats} \
        -O ${sampleID}.filtered.vcf.gz \
        --mitochondria-mode \
        ${extra_args}

    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp"  VariantFiltration \
        -V ${sampleID}.filtered.vcf.gz \
        -O ${sampleID}.mutect2.exclusionFiltered.${type}.vcf.gz \
        --apply-allele-specific-filters \
        --mask ${params.exclusion_sites} \
        --mask-name "exclusion_site"
    """
}
