process EMASE_BAM2EMASE {
    tag "$sampleID"

    cpus 1
    memory 200.GB
    time 6.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:v1.1.0-338c782'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.h5", mode: 'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.emase.h5"), emit: emase_h5

    script:
    """
    emase bam2emase -i ${bam} \
                -m ${params.transcripts_info} \
                -h ${params.gbrs_strain_list} \
                -o ${bam.baseName}.emase.h5
    """

    stub:
    """
    touch ${bam.baseName}.emase.h5
    """
}

/*
 Usage: emase bam2emase [OPTIONS]

 Convert BAM alignment files to EMASE format for allele-specific expression analysis

╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --alignment-files  -i      FILE     Input BAM file containing RNA-seq alignments, can separate files by "," or have multiple -i    │
│                                        [default: None]                                                                                │
│                                        [required]                                                                                     │
│ *  --haplotype-char   -h      TEXT     Haplotype identifiers (e.g., A,B,C,D). Can specify multiple times or comma-separated           │
│                                        [default: None]                                                                                │
│                                        [required]                                                                                     │
│ *  --locus-ids        -m      FILE     Transcript/locus information file [default: None] [required]                                   │
│    --output           -o      FILE     Output EMASE file (HDF5 format). Auto-generated if not specified [default: None]               │
│    --delim            -d      TEXT     Delimiter between transcript ID and haplotype in BAM file [default: _]                         │
│    --index-dtype              TEXT     Data type for matrix indices (advanced users only) [default: uint32]                           │
│    --data-dtype               TEXT     Data type for matrix values (advanced users only) [default: uint8]                             │
│    --verbose          -v      INTEGER  Increase verbosity (use multiple times for more detail) [default: 0]                           │
│    --help                              Show this message and exit.                                                                    │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
*/
