process HAPLOCHECK {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/haplocheck:1.3.3--h2a3209d_2'

    publishDir "${params.pubdir}/${sampleID + '/contamination'}", pattern: "*contamination.txt", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(tbi)

    output:
    tuple val(sampleID), file("*contamination.txt"), emit: contamination_txt
    tuple val(sampleID), file("*major_hg.txt"), emit: major_hg
    tuple val(sampleID), file("*minor_hg.txt"), emit: minor_hg
    tuple val(sampleID), file("*mean_het_major.txt"), emit: mean_het
    tuple val(sampleID), file("*mean_het_minor.txt"), emit: mean_het_minor
    
    script:
    """
    set -e
    PARENT_DIR="\$(dirname "${input_vcf}")"
    java -jar /usr/mtdnaserver/haplocheckCLI.jar "${PARENT_DIR}"

    sed 's/\"//g' output > output-noquotes

    grep "SampleID" output-noquotes > headers

    FORMAT_ERROR="Bad contamination file format"

    if [ `awk '{print \$2}' headers` != "Contamination" ]; then
      echo \$FORMAT_ERROR; exit 1
    fi
    if [ `awk '{print \$6}' headers` != "HgMajor" ]; then
      echo \$FORMAT_ERROR; exit 1
    fi
    if [ `awk '{print \$8}' headers` != "HgMinor" ]; then
      echo \$FORMAT_ERROR; exit 1
    fi
    if [ `awk '{print \$14}' headers` != "MeanHetLevelMajor" ]; then
      echo \$FORMAT_ERROR; exit 1
    fi
    if [ `awk '{print \$15}' headers` != "MeanHetLevelMinor" ]; then
      echo \$FORMAT_ERROR; exit 1
    fi

    grep -v "SampleID" output-noquotes > output-data
    awk -F "\t" '{print \$2}' output-data > ${sampleID}_contamination.txt
    awk -F "\t" '{print \$6}' output-data > ${sampleID}_major_hg.txt
    awk -F "\t" '{print \$8}' output-data > ${sampleID}_minor_hg.txt
    awk -F "\t" '{print \$14}' output-data > ${sampleID}_mean_het_major.txt
    awk -F "\t" '{print \$15}' output-data > ${sampleID}_mean_het_minor.txt
    """

}
