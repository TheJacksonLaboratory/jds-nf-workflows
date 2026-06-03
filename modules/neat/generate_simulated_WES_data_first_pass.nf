process GENERATE_SIMULATED_WES_DATA_FIRST_PASS {
    tag "$chunk_fasta.simpleName"


    cpus 2
    memory 50.GB
    time 8.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/neat:3.3'


    input:
    path chunk_fasta

    output:
    path('*vcf'),     emit: vcf
 

    script:
    prefix = params.gen_org=='mouse' ? "Mus_musculus.GRCm38" : "Homo_sapiens.GRCh38"
    """
    python /usr/local/bin/NEAT/gen_reads.py \
        -r ${chunk_fasta} \
        -R 150 \
        -o ${prefix}_simVar_60x_${chunk_fasta.simpleName} \
        -c 120 \
        -to 0.02 \
        -E 0.001 \
        --vcf \
        --pe 350 30 \
        -tr ${params.target_bed} \
        --no-fastq

    gunzip *vcf.gz
    """
}
