process OPTITYPE_RUN {
    tag "$sampleID"

    cpus 2
    memory 50.GB
    time 8.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/optitype:v1.5.1'

    publishDir "${params.pubdir}/hla_typing", pattern: '*.tsv', mode:'copy'  
    publishDir "${params.pubdir}/hla_typing", pattern: '*.pdf', mode:'copy'  

    input:
    tuple val(sampleID), path(fq_reads)

    output:
    path('*.tsv'), emit: tsv
    path('*.pdf'), emit: pdf, optional: true
 
    script:
    if (params.gen_org == 'human' && (params.workflow=='rnaseq' || params.workflow=='rna_fusion')){
      seq_type = "rna"
    } else {
      seq_type = "dna"
    } 
    paired_option = params.read_type == 'PE' ? "${fq_reads[0]} -i ${fq_reads[1]}" : "${fq_reads[0]}"
    """
    mkdir -p /flashscratch/${USER}/.config/matplotlib
    export MPLCONFIGDIR=/flashscratch/${USER}/.config/matplotlib
    optitype run -i ${paired_option} \
        --${seq_type} \
        -v \
        --outdir ./ \
	  -p ${sampleID}

    """
}
