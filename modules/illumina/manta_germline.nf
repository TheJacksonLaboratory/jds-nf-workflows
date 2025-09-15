process MANTA_CALL {
    tag "$sampleID"
    
    cpus = 8
    memory 80.GB
    time {bam.size() < 40.GB ? '10:00:00' : '24:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/manta:1.6.0--py27_0'

    input:
        tuple val(sampleID), path(bam), path(bai)
        tuple path(fasta), path(fai)
    
    output:
        tuple val(sampleID), path("${sampleID}_mantaDiploidSV.vcf"), emit: manta_sv

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    rm -rf mantaSVOut
    python /usr/local/bin/configManta.py \
        --runDir mantaSVOut \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        --callRegions ${params.callRegions}

    python ./mantaSVOut/runWorkflow.py \
    -m local -j ${task.cpus}

    mv mantaSVOut/results/variants/candidateSV.vcf.gz ./${sampleID}_mantaCandidate.vcf.gz
    mv mantaSVOut/results/variants/diploidSV.vcf.gz ./${sampleID}_mantaDiploidSV.vcf.gz
    gunzip ${sampleID}_mantaDiploidSV.vcf.gz
    """
}
