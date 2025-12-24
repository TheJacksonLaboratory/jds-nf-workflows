process STAR_ALIGN {
    tag "$sampleID"

    cpus 12
    memory 84.GB
    time 24.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*genome.sorted.ba*", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*transcript.sorted.ba*", mode:'copy'

    input:
    tuple val(sampleID), path(reads), val(read_length)
    val(rsem_ref_path)
    val(rsem_star_prefix)

    output:
    tuple val(sampleID), path("*.genome.sorted.bam"), path("*.genome.sorted.bam.bai"), emit: sorted_genomic_bam_bai
    tuple val(sampleID), path("*.transcript.sorted.bam"), path("*.transcript.sorted.bam.bai"), emit: sorted_transcript_bam_bai
    tuple val(sampleID), path("*final.out"), emit: star_log


    script:
    sort_command="samtools sort -@ 6 -m 5G -o ${sampleID}.transcript.sorted.bam ${sampleID}_Aligned.toTranscriptome.out.bam && samtools sort -@ 6 -m 5G -o ${sampleID}.genome.sorted.bam ${sampleID}_Aligned.out.bam"
    index_command="samtools index ${sampleID}.transcript.sorted.bam && samtools index ${sampleID}.genome.sorted.bam"

    read_length = read_length.toInteger()

    if( read_length >= 65 && read_length <= 85) {
        rsem_ref_files = file("${rsem_ref_path}/STAR/${rsem_star_prefix}_75/*").collect { "$it" }.join(' ')
    } else if( read_length >= 90 && read_length <= 110 ) {
        rsem_ref_files = file("${rsem_ref_path}/STAR/${rsem_star_prefix}_100/*").collect { "$it" }.join(' ')
    } else if( read_length >= 115 && read_length <= 135 ) {
        rsem_ref_files = file("${rsem_ref_path}/STAR/${rsem_star_prefix}_125/*").collect { "$it" }.join(' ')
    } else if( read_length >= 140 && read_length <= 160 ) {
        rsem_ref_files = file("${rsem_ref_path}/STAR/${rsem_star_prefix}_150/*").collect { "$it" }.join(' ')
    } else {
        log.info("\nUnsupported read length " + read_length + " in RSEM with STAR. RSEM will now fail gracefully.\n\n")
        rsem_ref_files = 'error'
    }

    """
    if [ "${rsem_ref_files}" = "error" ]; then exit 1; fi

    ln -s -f ${rsem_ref_files} .

    STAR \
    --genomeDir . \
    --readFilesIn ${reads}  \
    --runThreadN ${task.cpus} \
    --outFileNamePrefix ${sampleID}_ \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --genomeLoad NoSharedMemory \
    --outSAMtype BAM Unsorted \
    --quantMode TranscriptomeSAM
    
    ${sort_command}

    ${index_command}
    """
}

/*
NOTE: The STAR flag setting above are dictated by RSEM, which uses the ENCODE3 RNA-seq settings as a baseline

The codeblock from RSEM is shown below for reference:

    $command = "$star_path"."STAR" . 
                   ## ENCODE3 pipeline parameters
                   " --genomeDir $star_genome_path " .
                   " --outSAMunmapped Within " .
                   " --outFilterType BySJout " .
                   " --outSAMattributes NH HI AS NM MD " .
                   " --outFilterMultimapNmax 20 " .
                   " --outFilterMismatchNmax 999 " .
                   " --outFilterMismatchNoverLmax 0.04 " .
                   " --alignIntronMin 20 " .
                   " --alignIntronMax 1000000 " .
                   " --alignMatesGapMax 1000000 " .
                   " --alignSJoverhangMin 8 " .
                   " --alignSJDBoverhangMin 1 " .
                   " --sjdbScore 1 " .
                   " --runThreadN $nThreads " .
                   ##

                   ## different than ENCODE3 pipeline 
                   ## do not allow using shared memory
                   " --genomeLoad NoSharedMemory " .
                   ##

                   ## different than ENCODE3 pipeline, which sorts output BAM
                   ## no need to do it here to save time and memory 
                   " --outSAMtype BAM Unsorted " .
                   ##

                   ## unlike ENCODE3, we don"t output bedGraph files

                   " --quantMode TranscriptomeSAM ".
                   " --outSAMheaderHD \@HD VN:1.4 SO:unsorted ".

                   ## define output file prefix
                   " --outFileNamePrefix $imdName ";
                   ##

... ... ... 

    if ($mTime) { $time_start = time(); }

    &runCommand($command);

    if ($mTime) { $time_end = time(); $time_alignment = $time_end - $time_start; }

    $inpF = "$imdName.bam";

    if ( $star ) {
        my $star_tr_bam = $imdName . "Aligned.toTranscriptome.out.bam";
        rename $star_tr_bam, $inpF or die "can't rename $star_tr_bam to $inpF: $!\n";
        rmdir $imdName . "_STARtmp/";
        my $star_genome_bam = $imdName . "Aligned.out.bam";
        my $rsem_star_genome_bam = $sampleName.'.STAR.genome.bam';
        if ( $star_output_genome_bam ) {
            rename $star_genome_bam, $rsem_star_genome_bam or die "can't move $star_genome_bam to $rsem_star_genome_bam: $!\n";
        } else {
            unlink $star_genome_bam or die "can't remove $star_genome_bam: $!\n";
        }
        rename $imdName."Log.final.out", $sampleName.".log" or die "Cannot rename ${imdName}Log.final.out to $sampleName.log: $!\n";
    }
*/

