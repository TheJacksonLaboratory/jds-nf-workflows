process HAPLOCHECK {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/haplocheck:1.3.3--h2a3209d_2'

    publishDir "${params.pubdir}/${sampleID + '/mt_callers/contamination'}", pattern: "*haplocheck_output.txt", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(tbi)

    output:
    tuple val(sampleID), path("*haplocheck_output.txt"), emit: contam_report
    tuple val(sampleID), path("*contamination.txt"), emit: contam_value
    tuple val(sampleID), env("contam_value"), emit: contam
    
    script:
    """
    set -e
    haplocheck --raw --out=output ${vcf} 

    sed 's/\"//g' output.raw.txt > ${sampleID}.haplocheck_output.txt

    # Extract the second row (exclude header) and check if \$2 == 'YES'
    contam_value=\$(awk 'NR==2 && \$2=="YES" {print \$3}' ${sampleID}.haplocheck_output.txt)
    if [ -n "\$contam_value" ]; then
      echo "\$contam_value" > ${sampleID}.contamination.txt
    else
      echo "0" > ${sampleID}.contamination.txt
    fi 
    """

}

/*
In version 1.3.3 of Haplocheck, the 'Overall Contamination' is set by: 

	private double calcOverallLevel(double major, double minor) {

		if (major > 0) {
			return 1 - major;
		} else {
			return minor;
		}

	}

Therefore, we can use that value, and omit the following section of
code that was provided in: https://github.com/broadinstitute/gatk/blob/4.1.8.0/scripts/mitochondria_m2_wdl/AlignAndCall.wdl

...
set -e
PARENT_DIR="$(dirname "~{input_vcf}")"
java -jar /usr/mtdnaserver/haplocheckCLI.jar "${PARENT_DIR}"

sed 's/\"//g' output > output-noquotes

grep "SampleID" output-noquotes > headers
FORMAT_ERROR="Bad contamination file format"
if [ `awk '{print $2}' headers` != "Contamination" ]; then
  echo $FORMAT_ERROR; exit 1
fi
if [ `awk '{print $6}' headers` != "HgMajor" ]; then
  echo $FORMAT_ERROR; exit 1
fi
if [ `awk '{print $8}' headers` != "HgMinor" ]; then
  echo $FORMAT_ERROR; exit 1
fi
if [ `awk '{print $14}' headers` != "MeanHetLevelMajor" ]; then
  echo $FORMAT_ERROR; exit 1
fi
if [ `awk '{print $15}' headers` != "MeanHetLevelMinor" ]; then
  echo $FORMAT_ERROR; exit 1
fi

grep -v "SampleID" output-noquotes > output-data
awk '{print $2}' output-data > contamination.txt
awk '{print $6}' output-data > major_hg.txt
awk '{print $8}' output-data > minor_hg.txt
awk '{print $14}' output-data > mean_het_major.txt
awk '{print $15}' output-data > mean_het_minor.txt

...

output {
File contamination_file = "output-noquotes"
String hasContamination = read_string("contamination.txt") 
String major_hg = read_string("major_hg.txt")
String minor_hg = read_string("minor_hg.txt")
Float major_level = read_float("mean_het_major.txt")
Float minor_level = read_float("mean_het_minor.txt")
}

... 

Float hc_contamination = if run_contamination && hasContamination == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
Float max_contamination = if defined(verifyBamID) && verifyBamID > hc_contamination then verifyBamID else hc_contamination

*/
