process SNPEFF {
    tag "$sampleID"

    cpus 1
    memory  16.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d'

    // chip_somatic counterpart of snpeff_snpsift/snpeff_snpeff.nf (unchanged),
    // trimmed to the one case chip_somatic actually uses (indel_snp == 'SNP' or
    // 'INDEL', both species -- see config/chip_somatic.config for gen_ver/
    // snpEff_config resolved per species). No publishDir: intermediate only.

    input:
    tuple val(sampleID), file(vcf)
    val(indel_snp)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    output_suffix = indel_snp == 'INDEL' ? 'INDEL_snpeff.vcf' : 'SNP_snpeff.vcf'
    """
    java -Djava.io.tmpdir=./ -Xmx${my_mem}G -jar /opt/snpEff/snpEff.jar \
    ${params.gen_ver} \
    -c ${params.snpEff_config} \
    -o vcf \
    -noStats \
    ${vcf} > ${sampleID}_${output_suffix}
    """
}
