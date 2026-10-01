process GATHER_SUMMARIES {
    tag "all samples"

    input:
    path summary_files

    output:
    path "all_summaries.txt", emit: gathered_file

    script:
    """
    cat ${summary_files} > all_summaries.txt
    """
}