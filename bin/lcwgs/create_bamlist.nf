def create_bamlist( bams, downsample_value ) {
    bams.collectFile(name = "bamlist.txt", newline: true).set { bam_list }
    return [bam_list, downsample_value]
}