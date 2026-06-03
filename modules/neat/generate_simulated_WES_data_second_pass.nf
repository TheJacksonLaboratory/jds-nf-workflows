process GENERATE_SIMULATED_WES_DATA_SECOND_PASS {
    tag "$chunk_fasta.simpleName"


    cpus 2
    memory 50.GB
    time 8.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/neat:3.3'

    publishDir "${params.pubdir}/results/simulated_60x_individual_chr", pattern: '*fq.gz', mode:'copy', enabled: params.workflow=='generate_wes_simreads'  


    input:
    path chunk_fasta
    path int_vcf

    output:
    path('*1.fq.gz'), emit: fq1, optional: true
    path('*2.fq.gz'), emit: fq2, optional: true
    path('*vcf'),     emit: vcf
 

    script:
    prefix = params.gen_org=='mouse' ? "Mus_musculus.GRCm38" : "Homo_sapiens.GRCh38"
    target_bed = params.gen_org=='mouse' ? params.target_bed : params.target_sort_bed

    """
    python /usr/local/bin/NEAT/gen_reads.py \
        -r ${chunk_fasta} \
        -R 150 \
        -o ${prefix}_simVar_60x_${chunk_fasta.simpleName} \
        -c 120 \
        -M 0 \
        -to 0.02 \
        -E 0.001 \
        --vcf \
        -v ${int_vcf} \
        --pe 350 30 \
        -tr ${target_bed} 

    gunzip *vcf.gz
        
    """
}
