process PORYGON {
    tag "used python. It was super effective."

    cpus 1
    time '00:05:00'
    memory  1.GB

    script:
    """
    python ${moduleDir}/bin/porygon.py
    """
}

// This module is for testing module / nextflow failures only. It is not intended for production level pipeline use. 
