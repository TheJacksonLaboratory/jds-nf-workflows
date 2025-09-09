process FASTP {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "${sampleID}.fastp*.{json,html,log}", mode:'copy'


    input:
    tuple val(sampleID), path(reads)
    path  adapter_fasta
    val   save_trimmed_fail
    val   save_merged

    output:
    tuple val(sampleID), path('*.fastp.fastq.gz') , optional:true, emit: reads
    tuple val(sampleID), path('*.json')           , emit: json
    tuple val(sampleID), path('*.html')           , emit: html
    tuple val(sampleID), path('*.log')            , emit: log
    tuple val(sampleID), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(sampleID), path('*.merged.fastq.gz'), optional:true, emit: reads_merged
 

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sampleID}"
    def adapter_list = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
    def fail_fastq = save_trimmed_fail && params.read_type == 'SE' ? "--failed_out ${prefix}.fail.fastq.gz" : save_trimmed_fail && !params.read_type == 'SE' ? "--unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    // Use single ended for interleaved. Add --interleaved_in in config.
    if ( task.ext.args?.contains('--interleaved_in') ) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        fastp \\
            --stdout \\
            --in1 ${prefix}.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $args \\
            2> ${prefix}.fastp.log \\
        | gzip -c > ${prefix}.fastp.fastq.gz

        """
    } else if (params.read_type == 'SE') {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        fastp \\
            --in1 ${prefix}.fastq.gz \\
            --out1  ${prefix}.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $args \\
            2> ${prefix}.fastp.log

        """
    } else {
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.fastp.fastq.gz \\
            --out2 ${prefix}_2.fastp.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2> ${prefix}.fastp.log

        """
    }
}
