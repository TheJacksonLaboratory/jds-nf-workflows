process PAV {

    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "24:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'becklab/pav:2.4.6'

    publishDir "${params.pubdir}/${sampleID + '/alignments'}", pattern: "${sampleID}.pbmm2.aligned.bam*", mode: "copy"

    input:
        tuple val(sampleID), path(fq)
        path(ref_fa)

    output:
        val(sampleID), emit: pav_output
    
    script:
    
    """
    echo "Making config.json for PAV"
    echo '{ "reference": "${ref_fa}" }' > config.json

    echo "Making assemblies.tsv for PAV"
    echo -e "NAME\\tHAP_unphased\\n${sampleID}\\t${fq}" > assemblies.tsv

    echo "Running PAV"
    export HOME=$PWD
    export TMPDIR=$PWD
    /opt/pav/files/docker/run
    """
}
