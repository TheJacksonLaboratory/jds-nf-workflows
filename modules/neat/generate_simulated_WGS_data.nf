process GENERATE_SIMULATED_WGS_DATA {
    tag "$chunk_fasta.simpleName"

    cpus 2
    memory 50.GB
    time 8.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/neat:3.3'

    publishDir "${params.pubdir}/reads/individual_chr", pattern: '*fq.gz', mode:'copy'  
    publishDir "${params.pubdir}/gold_truth_vcf", pattern: '*golden.vcf', mode:'copy'  

    input:
    path chunk_fasta

    output:
    path('*1.fq.gz'), emit: fq1
    path('*2.fq.gz'), emit: fq2, optional: true
    path('*vcf'),     emit: vcf
 
    script:
    prefix = params.gen_org=='mouse' ? "Mus_musculus.GRCm38" : "Homo_sapiens.GRCh38"
    paired_end = params.read_type=='PE' ? "--pe ${params.insert_size} ${params.insert_size_sd}" : ""
    """
    python /usr/local/bin/NEAT/gen_reads.py \
        -r ${chunk_fasta} \
        -R ${params.read_length} \
        -o ${prefix}_simVar_10x_${chunk_fasta.simpleName} \
	    -c ${params.coverage} \
        -E ${params.error_rate} \
        ${paired_end} \
        --vcf

    gunzip *vcf.gz
    """
}
