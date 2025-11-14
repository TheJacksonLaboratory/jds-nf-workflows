process MIRTRACE_RUN {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '08:00:00'


    container 'quay.io/biocontainers/mirtrace:1.0.1--hdfd78af_1'

    publishDir "${params.pubdir}", pattern: "mirtrace/*", mode:'copy'

    input:
    tuple val(adapter), val(ids), path(reads)

    output:
    path "mirtrace/*"  , emit: mirtrace

    when:
    task.ext.when == null || task.ext.when

    script:
    // mirtrace protocol defaults to 'params.protocol' if not set
    def primer = adapter ? "--adapter ${adapter}" : ""
    def protocol = params.protocol == 'custom' ? '' : "--protocol $params.protocol"
    def java_mem = ''
    if(task.memory){
        tmem = task.memory.toBytes()
        java_mem = "-Xms${tmem} -Xmx${tmem}"
    }
    def config_lines = [ids,reads]
    .transpose()
    .collect({ id, path -> "echo '${path},${id}' >> mirtrace_config" })
    """
    export mirtracejar=\$(dirname \$(which mirtrace))

    ${config_lines.join("\n    ")}

    java $java_mem -jar \$mirtracejar/mirtrace.jar --mirtrace-wrapper-name mirtrace qc  \\
        --species $params.mirtrace_species \\
        $primer \\
        $protocol \\
        --config mirtrace_config \\
        --write-fasta \\
        --output-dir mirtrace \\
        --force

    """

}
