process PBMM2_CALL {

    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "24:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern:"*.bam", mode:'copy', enabled: (params.merge_inds == false) ? true : false

    input:
        tuple val(sampleID), path(fq)
        path(pbmm2_index)

    output:
        tuple val(sampleID), file("${sampleID}.pbmm2.aligned.bam"), file("${sampleID}.pbmm2.aligned.bam.bai"), emit: pbmm2_bam 
    
    script:
        if (params.pbmode == "CCS")
            """
            pbmm2 align ${pbmm2_index} ${fq} ${sampleID}.pbmm2.aligned.bam --preset CCS --sort -j ${task.cpus}
            """
        else if (params.pbmode == "CLR")
            """
            pbmm2 align ${pbmm2_index} ${fq} ${sampleID}.pbmm2.aligned.bam --median-filter --sort -j ${task.cpus}
            """
}
