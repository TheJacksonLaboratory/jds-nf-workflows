process MIRTOP_QUANT {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:63be019f50581cc5dfe4fc0f73ae50f2d4d661f7-0' 

    publishDir "${params.pubdir}", pattern: "mirtop/mirtop*", mode:'copy'
    publishDir "${params.pubdir}", pattern: "mirtop/stats/*", mode:'copy'


    input:
    path ("bams/*")
    path hairpin
    path gtf

    output:
    path "mirtop/mirtop.gff"        , emit: mirtop_gff
    path "mirtop/mirtop.tsv"        , emit: mirtop_table
    path "mirtop/mirtop_rawData.tsv", emit: mirtop_rawdata
    path "mirtop/stats/*"           , emit: logs


    script:
    def filter_species = params.mirgenedb ? params.mirgenedb_species : params.mirtrace_species
    """
    #Cleanup the GTF if mirbase html form is broken
    GTF="$gtf"
    sed 's/&gt;/>/g' \$GTF | sed 's#<br>#\\n#g' | sed 's#</p>##g' | sed 's#<p>##g' | sed -e :a -e '/^\\n*\$/{\$d;N;};/\\n\$/ba' > \${GTF}_html_cleaned.gtf
    mirtop gff --hairpin $hairpin --gtf \${GTF}_html_cleaned.gtf -o mirtop --sps $filter_species ./bams/*
    mirtop counts --hairpin $hairpin --gtf \${GTF}_html_cleaned.gtf -o mirtop --sps $filter_species --add-extra --gff mirtop/mirtop.gff
    mirtop export --format isomir --hairpin $hairpin --gtf \${GTF}_html_cleaned.gtf --sps $filter_species -o mirtop mirtop/mirtop.gff
    mirtop stats mirtop/mirtop.gff --out mirtop/stats
    mv mirtop/stats/mirtop_stats.log mirtop/stats/full_mirtop_stats.log

    """

}
