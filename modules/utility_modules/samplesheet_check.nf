process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    container 'quay.io/biocontainers/python:3.8.3'

    input:
    path(samplesheet)

    output:
    path('*.csv'), emit: csv

    script:
    """
    ${projectDir}/bin/smrnaseq/check_samplesheet.py \
    $samplesheet \
    samplesheet.valid.csv
    """
}

