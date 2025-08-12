process PICARD_LIFTOVERVCF_MERGEVCF {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
    
    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*.vcf*", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(vcf), path(shifted_vcf)

    output:
    tuple val(sampleID), file("*final.vcf"), emit: vcf
    tuple val(sampleID), file("*final.vcf.idx"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp LiftoverVcf \
      I=${shifted_vcf} \
      O=${shifted_vcf.baseName}.shifted_back.vcf \
      R=${params.mt_fasta} \
      CHAIN=${params.shift_back_chain} \
      REJECT=${shifted_vcf.baseName}.rejected.vcf

    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp MergeVcfs \
      I=${shifted_vcf.baseName}.shifted_back.vcf \
      I=${vcf} \
      O=${shifted_vcf.baseName}.final.vcf
    """
}
