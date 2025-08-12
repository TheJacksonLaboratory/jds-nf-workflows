process PICARD_MT_COVERAGEATEVERYBASE {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
    
    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*per_base_coverage.tsv", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(bam), path(bai), path(shifted_bam), path(shifted_bai)

    output:
    tuple val(sampleID), file("*per_base_coverage.tsv"), emit: coverage

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    set -e

    # Run Picard CollectHsMetrics for non-control region
    mkdir -p tmp
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp CollectHsMetrics \
      I=${bam} \
      R=${params.mt_fasta} \
      PER_BASE_COVERAGE=non_control_region.tsv \
      O=non_control_region.metrics \
      TI=${params.non_control_region_interval_list} \
      BI=${params.non_control_region_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    # Run Picard CollectHsMetrics for control region (shifted)
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp CollectHsMetrics \
      I=${shifted_bam} \
      R=${params.mt_shifted_fasta} \
      PER_BASE_COVERAGE=control_region_shifted.tsv \
      O=control_region_shifted.metrics \
      TI=${params.control_region_shifted_interval_list} \
      BI=${params.control_region_shifted_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    # Write the R script to a file
    cat > merge_coverage.R <<'EOF'
    shift_back <- function(x) {
      if (x < 8570) {
        return(x + 8000)
      } else {
        return(x - 8569)
      }
    }

    control_region_shifted <- read.table("control_region_shifted.tsv", header=TRUE, sep="\\t")
    control_region_shifted\$pos <- sapply(control_region_shifted\$pos, shift_back)

    beginning <- subset(control_region_shifted, pos < 8000)
    end <- subset(control_region_shifted, pos > 8000)

    non_control_region <- read.table("non_control_region.tsv", header=TRUE, sep="\\t")
    combined_table <- rbind(beginning, non_control_region, end)
    write.table(combined_table, paste0("${sampleID}_per_base_coverage.tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\\t")
    EOF

    # Run the R script
    Rscript merge_coverage.R
    """
}
