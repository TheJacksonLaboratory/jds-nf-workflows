process MAKE_PERM_VALUE_FILE {

    tag "$phenotype"
    
    cpus 1
    memory '5 GB'
    time 1.hour
    
    container 'ubuntu:20.04'

    input:
    tuple val(id), val(phenotype)

    output:
    tuple val(id), val(phenotype), path('*_scan1perms.txt'), emit: perm_files

    script:

    """
    echo "perm" > ${phenotype}_scan1perms.txt
    echo ${params.perm_threshold} >> ${phenotype}_scan1perms.txt
    """

}