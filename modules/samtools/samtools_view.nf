process SAMTOOLS_VIEW {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools_view' }", pattern:"*.bam", mode:'copy', enabled: params.keep_intermediate

  input:
      tuple val(sampleID), file(sam)
      val(view_string)
      val(filename)

  output:
      tuple val(sampleID), file("*.bam"), emit: bam

  script:
    if (params.workflow == 'rna_fusion')
    """
    awk '
      BEGIN { FS="\\t"; OFS="\\t" }
      /^@/ { print; next }
      {
        # Exclude records with POS (\$4) or PNEXT (\$8) beyond BAM coordinate limits
        if (\$4 > 536870911 || \$8 > 536870911) next

        # Sanitize integer overflow in TLEN (\$9) to 0 (valid SAM spec for unknown insert size)
        if (\$9 > 2147483647 || \$9 < -2147483647) {
          \$9 = 0
        }

        print
      }
    ' ${sam} | samtools view ${view_string} - > ${sampleID}_${filename}.bam
    """
    else
    """
    samtools view ${view_string} ${sam} > ${sampleID}_${filename}.bam
    """
  
  stub:
    """
    ${sampleID}_${filename}.bam
    """
}

/* Regarding the rna_fusion workflow:

The input sam file from STAR_SQUID contains chimeric reads.
In cases where STAR aligns complex chimeric or distant split reads, its internal math calculates an infinite or wrapped-around distance, outputting 2^{64}-1 to TLEN.
The value 2^{64} - 1 = 18446744073709551615, which is outside the valid range for a signed 32-bit integer.
The TLEN field is defined as a signed 32-bit integer in the SAM specification, which means that the maximum valid value for TLEN is 2^31 - 1 = 2147483647. 
When samtools view encounters a TLEN value outside this range, it will cause an integer overflow error.
To prevent this, we sanitize the TLEN field by setting it to 0 for any values outside the valid range. This ensures that samtools view can process the SAM file without encountering integer overflow errors.

We further filter out any records with POS or PNEXT values greater than 536870911, which is the maximum valid coordinate for BAM files. This prevents any potential issues with invalid coordinates in the output BAM file.
*/