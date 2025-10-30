process PICARD_MERGEBAMALIGNMENT {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    input:
    tuple val(sampleID), path(mapped_bam), path(unmapped_bam)
    val(type)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    reference = type == 'shifted_mt' ? params.mt_shifted_fasta : params.mt_fasta

    """
    mkdir -p tmp
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp MergeBamAlignment \
    VALIDATION_STRINGENCY=SILENT \
    EXPECTED_ORIENTATIONS=FR \
    ATTRIBUTES_TO_RETAIN=X0 \
    ATTRIBUTES_TO_REMOVE=NM \
    ATTRIBUTES_TO_REMOVE=MD \
    ALIGNED_BAM=${mapped_bam} \
    UNMAPPED_BAM=${unmapped_bam} \
    OUTPUT=${mapped_bam.baseName}.mba.bam \
    REFERENCE_SEQUENCE=${reference} \
    PAIRED_RUN=true \
    SORT_ORDER="unsorted" \
    IS_BISULFITE_SEQUENCE=false \
    ALIGNED_READS_ONLY=false \
    CLIP_ADAPTERS=false \
    MAX_RECORDS_IN_RAM=2000000 \
    ADD_MATE_CIGAR=true \
    MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    PROGRAM_RECORD_ID="bwamem" \
    PROGRAM_GROUP_VERSION="0.7.17" \
    PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 2 -Y" \
    PROGRAM_GROUP_NAME="bwamem" \
    UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
    ALIGNER_PROPER_PAIR_FLAGS=true \
    UNMAP_CONTAMINANT_READS=true \
    ADD_PG_TAG_TO_READS=false
    """
}
