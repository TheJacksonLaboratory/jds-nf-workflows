process RSEM_ALIGNMENT_EXPRESSION {
    tag "$sampleID"

    cpus 12
    memory 90.GB
    time 24.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0'

    publishDir path: {"${params.pubdir}/${sampleID + '/stats'}"}, pattern: "*stats", mode:'copy', enabled: params.rsem_aligner == "bowtie2"
    publishDir path: {"${params.pubdir}/${sampleID + '/stats'}"}, pattern: "*final.out", mode:'copy', enabled: params.rsem_aligner == "star"
    publishDir path: {"${params.pubdir}/${sampleID}"}, pattern: "*results*", mode:'copy'
    publishDir path: {"${params.pubdir}/${sampleID + '/bam'}"}, pattern: "*genome.sorted.ba*", mode:'copy'
    publishDir path: {"${params.pubdir}/${sampleID + '/bam'}"}, pattern: "*transcript.sorted.ba*", mode:'copy'

    input:
    tuple val(sampleID), path(reads), val(strand_setting), val(read_length)
    val(rsem_ref_path)
    val(rsem_star_prefix)
    val(rsem_ref_prefix)

    output:
    path "*stats"
    path "*results*"
    tuple val(sampleID), path("rsem_aln_*.stats"), emit: rsem_stats
    tuple val(sampleID), path("*.stat/*.cnt"), emit: rsem_cnt
    tuple val(sampleID), path("*genes.results"), emit: rsem_genes
    tuple val(sampleID), path("*isoforms.results"), emit: rsem_isoforms
    tuple val(sampleID), path("*.genome.bam"), emit: bam
    tuple val(sampleID), path("*.transcript.bam"), emit: transcript_bam
    tuple val(sampleID), path("*.genome.sorted.bam"), emit: sorted_genomic_bam
    tuple val(sampleID), path("*.genome.sorted.bam.bai"), emit: sorted_genomic_bam_bai, optional: true
    tuple val(sampleID), path("*.transcript.sorted.bam"), emit: sorted_transcript_bam
    tuple val(sampleID), path("*.transcript.sorted.bam.bai"), emit: sorted_transcript_bam_bai, optional: true
    tuple val(sampleID), path("*final.out"), emit: star_log, optional: true

    /* 
        Note: For large genomes, the index step may fail. 
              skip_index may be used to disable index creation. 
              This is why sorted_genomic_bam and sorted_transcript_bam are marked as optional outputs as they may not be created if the index step isn't run.
    */
    
    script:

    def prob = ''
    if (strand_setting == "reverse_stranded") {
        prob="--forward-prob 0"
    }

    if (strand_setting == "forward_stranded") {
        prob="--forward-prob 1"
    }

    if (strand_setting == "non_stranded") {
        prob="--forward-prob 0.5"
    }

    def frag = ''
    def stype = ''
    def trimmedfq = ''
    if (params.read_type == "PE"){
        stype="--paired-end"
        trimmedfq="${reads[0]} ${reads[1]}"
    }
    if (params.read_type == "SE"){
        frag="--fragment-length-mean ${params.fragment_length_mean} --fragment-length-sd ${params.fragment_length_sd}"
        trimmedfq="${reads[0]}"
    }

    def rsem_ref_files
    def outbam
    def seed_length
    def sort_command
    def index_command
    def intermediate
    def star_log
    def readFilesCommand

    if (params.rsem_aligner == "bowtie2"){
        
        rsem_ref_files = files("${rsem_ref_path}/bowtie2/*").collect { it -> "$it" }.join(' ')

        outbam="--output-genome-bam"
        seed_length="--seed-length ${params.seed_length}"
        sort_command="samtools sort -@ 6 -m 5G -o ${sampleID}.transcript.sorted.bam ${sampleID}.transcript.bam && samtools sort -@ 6 -m 5G -o ${sampleID}.genome.sorted.bam ${sampleID}.genome.bam"
        index_command="samtools index ${sampleID}.transcript.sorted.bam && samtools index ${sampleID}.genome.sorted.bam"
        intermediate=''
        star_log=''
        readFilesCommand=''
    }

    if (params.rsem_aligner == "star") {
        outbam="--star-output-genome-bam"
        seed_length=""
        sort_command="samtools sort -@ 6 -m 5G -o ${sampleID}.transcript.sorted.bam ${sampleID}.transcript.bam && samtools sort -@ 6 -m 5G -o ${sampleID}.STAR.genome.sorted.bam ${sampleID}.STAR.genome.bam"
        index_command="samtools index ${sampleID}.transcript.sorted.bam && samtools index ${sampleID}.STAR.genome.sorted.bam"
        intermediate='--keep-intermediate-files'
        star_log="cp ${sampleID}/${sampleID}.temp/*.final.out ${sampleID}/stats/${sampleID}.STAR.Log.final.out && rm -r ${sampleID}/${sampleID}.temp"

        def read_len = read_length.toInteger()

        if( read_len >= 45 && read_len <= 60) {
            rsem_ref_files = files("${rsem_ref_path}/STAR/${rsem_star_prefix}_50/*").collect { it -> "$it" }.join(' ')
        } else if( read_len >= 65 && read_len <= 85) {
            rsem_ref_files = files("${rsem_ref_path}/STAR/${rsem_star_prefix}_75/*").collect { it -> "$it" }.join(' ')
        } else if( read_len >= 90 && read_len <= 110 ) {
            rsem_ref_files = files("${rsem_ref_path}/STAR/${rsem_star_prefix}_100/*").collect { it -> "$it" }.join(' ')
        } else if( read_len >= 115 && read_len <= 135 ) {
            rsem_ref_files = files("${rsem_ref_path}/STAR/${rsem_star_prefix}_125/*").collect { it -> "$it" }.join(' ')
        } else if( read_len >= 140 && read_len <= 160 ) {
            rsem_ref_files = files("${rsem_ref_path}/STAR/${rsem_star_prefix}_150/*").collect { it -> "$it" }.join(' ')
        } else {
            log.info("\nUnsupported read length " + read_len + " in RSEM with STAR. RSEM will now fail gracefully.\n\n")
            rsem_ref_files = 'error'
        }

        readFilesCommand = reads[0].toString().endsWith('.gz') ? '--star-gzipped-read-file' : ''

    }

    index_command = params.skip_index ? '' : index_command
    
    /* 
        Note: For large genomes, the index step may fail. 
              skip_index may be used to disable index creation. 
              Sorting is still valid.
    */
    
    """
    if [ "${rsem_ref_files}" = "error" ]; then exit 1; fi

    ln -s -f ${rsem_ref_files} . 

    rsem-calculate-expression -p $task.cpus \
    ${prob} \
    ${stype} \
    ${frag} \
    --${params.rsem_aligner} \
    ${readFilesCommand} \
    --append-names \
    ${seed_length} \
    ${outbam} \
    ${trimmedfq} \
    ${rsem_ref_prefix} \
    ${sampleID} \
    ${intermediate} \
    2> rsem_aln_${sampleID}.stats

    ${star_log}

    ${sort_command}

    ${index_command}
    """
}
