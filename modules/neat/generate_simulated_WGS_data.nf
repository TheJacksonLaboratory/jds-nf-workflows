process GENERATE_SIMULATED_WGS_DATA {
    tag "$chunk_fasta.simpleName"


    cpus 2
    memory 50.GB
    time 8.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/neat:3.3'

    publishDir "${params.pubdir}/results/simulated_10x_individual_chr", pattern: '*fq.gz', mode:'copy'  
    publishDir "${params.pubdir}/results/gold_truth_vcf", pattern: '*golden.vcf', mode:'copy'  


    input:
    path chunk_fasta

    output:
    path('*1.fq.gz'), emit: fq1
    path('*2.fq.gz'), emit: fq2
    path('*vcf'),     emit: vcf
 

    script:
    """
    python /usr/local/bin/NEAT/gen_reads.py \
        -r ${chunk_fasta} \
        -R 150 \
        -o Mus_musculus.GRCm38_simVar_10x_${chunk_fasta.simpleName} \
	-c 10 \
        -E 0.001 \
        --vcf \
        --pe 350 30
    
    gunzip *vcf.gz
    """

}
