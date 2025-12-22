process UMITOOLS_PREPAREFORRSEM {
    tag "$sampleID"

    cpus 1
    memory 64.GB
    time 1.h

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/biocontainers/umi_tools:1.1.6--py311haab0aaa_0"

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.rsem_prep.bam", mode:'copy'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path('*.bam'), emit: bam
    tuple val(sampleID), path('*.log'), emit: log

    script:
    """
    umi_tools prepare-for-rsem \
        --stdin=${bam} \
        --stdout=${bam.baseName}.rsem_prep.bam \
        --log=${bam.baseName}.prepare_for_rsem.log \
    """
}

/*
prepare_for_rsem - make output from dedup or group compatible with RSEM

Usage: umi_tools prepare_for_rsem [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr

For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

Options:
  --version             show program's version number and exit

  RSEM preparation specific options:
    --tags=TAGS         Comma-seperated list of tags to transfer from read1 to
                        read2
    --sam               input and output SAM rather than BAM

  input/output options:
    -I FILE, --stdin=FILE
                        file to read stdin from [default = stdin].
    -L FILE, --log=FILE
                        file with logging information [default = stdout].
    -E FILE, --error=FILE
                        file with error information [default = stderr].
    -S FILE, --stdout=FILE
                        file where output is to go [default = stdout].
    --temp-dir=FILE     Directory for temporary files. If not set, the bash
                        environmental variable TMPDIR is used[default = None].
    --log2stderr        send logging information to stderr [default = False].
    --compresslevel=COMPRESSLEVEL
                        Level of Gzip compression to use. Default (6)
                        matchesGNU gzip rather than python gzip default (which
                        is 9)

  profiling options:
    --timeit=TIMEIT_FILE
                        store timeing information in file [none].
    --timeit-name=TIMEIT_NAME
                        name in timing file for this class of jobs [all].
    --timeit-header     add header for timing information [none].

  common options:
    -v LOGLEVEL, --verbose=LOGLEVEL
                        loglevel [1]. The higher, the more output.
    -h, --help          output short help (command line options only).
    --help-extended     Output full documentation
    --random-seed=RANDOM_SEED
                        random seed to initialize number generator with
                        [none].
*/