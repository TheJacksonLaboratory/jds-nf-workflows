process GENERATE_SIMULATED_WGS_DATA {
    tag "$chunk_fasta.simpleName"

    cpus 2
    memory 50.GB
    time 8.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/neat:3.3'

    publishDir "${params.pubdir}/reads/individual_chr", pattern: '*fq.gz', mode:'copy'  
    publishDir "${params.pubdir}/gold_truth_vcf/individual_chr", pattern: '*golden.vcf', mode:'copy'  

    input:
    path chunk_fasta

    output:
    path('*1.fq.gz'), emit: fq1
    path('*2.fq.gz'), emit: fq2, optional: true
    path('*vcf'),     emit: vcf
 
    script:
    paired_end = params.read_type=='PE' ? "--pe ${params.insert_size} ${params.insert_size_sd}" : ""
    mutation_rate = params.mutation_rate != null ? "-M ${params.mutation_rate}" : "" // NEAT scales mutation rate, if not user set use internal tool defaults.
    """
    python /usr/local/bin/NEAT/gen_reads.py \
        -r ${chunk_fasta} \
        -R ${params.read_length} \
        -o ${params.sampleID}_simVar_${params.coverage}x_${chunk_fasta.simpleName} \
	    -c ${params.coverage} \
        -E ${params.error_rate} \
        ${paired_end} \
        ${mutation_rate} \
        --vcf

    gunzip *vcf.gz
    """
}
