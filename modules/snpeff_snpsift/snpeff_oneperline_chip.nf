process SNPEFF_ONEPERLINE {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '00:10:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    // chip_somatic counterpart of snpeff_snpsift/snpeff_oneperline.nf (unchanged).
    // Same script (${moduleDir}/bin/vcfEffOnePerLine.pl resolves to the SAME
    // shared bin/ next to this file, since this module lives in the same
    // modules/snpeff_snpsift/ directory -- no need to duplicate the perl
    // script). No container, matching the original: runs on the node's own
    // perl, same as WES already does.

    input:
    tuple val(sampleID), file(vcf)
    val(indel_snp)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    output_suffix = indel_snp == 'INDEL' ? 'INDEL_snpeff.vcf' : 'SNP_snpeff.vcf'
    """
    cat ${vcf} | perl ${moduleDir}/bin/vcfEffOnePerLine.pl > ${sampleID}_oneperline_${output_suffix}
    """
}
