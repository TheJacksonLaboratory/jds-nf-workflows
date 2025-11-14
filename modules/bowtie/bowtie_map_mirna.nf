process BOWTIE_MAP_SEQ {
    tag "$sampleID"

    cpus 8
    memory 30.GB
    time 15.hour


    container    'quay.io/biocontainers/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:40128b496751b037e2bd85f6789e83d4ff8a4837-0'

    input:
    tuple val(sampleID), path(reads)
    path index

    output:
    tuple val(sampleID), path("*bam")           , emit: bam
    tuple val(sampleID), path('unmapped/*fq.gz'), emit: unmapped


    script:
    cmd = reads =~ /gz/ ?  "zcat" : "cat"
    """
    INDEX=`find -L ./ -name "*.3.ebwt" | sed 's/.3.ebwt//'`
    bowtie \\
         -x \$INDEX \\
        -q <(${cmd} $reads) \\
        -p ${task.cpus} \\
        -t \\
        -k 50 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        --un ${sampleID}_unmapped.fq -S > ${sampleID}.sam

    samtools view -bS ${sampleID}.sam > ${sampleID}.bam

    if [ ! -f  "${sampleID}_unmapped.fq" ]
    then
        touch ${sampleID}_unmapped.fq
    fi
    gzip ${sampleID}_unmapped.fq
    mkdir unmapped
    mv  ${sampleID}_unmapped.fq.gz  unmapped/.

    """

}
