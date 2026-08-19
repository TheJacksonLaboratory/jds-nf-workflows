process CHECK_STRANDEDNESS {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time '1:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/how-are-we-stranded-here:v1.0.1-e6ce74d'

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern:"*_strandedness.txt", mode:'copy'

    input:
    tuple val(sampleID), path(reads)
    path(gtf)

    output:
    tuple val(sampleID), env("STRAND"), emit: strand_setting
    tuple val(sampleID), path('*_strandedness.txt'), emit: strandedness_report

    script:
    paired = params.read_type == 'PE' ? "-r2 ${reads[1]}" : ''

    // In cases where the GTF is provided as a gzipped file, we need to unzip it for use in the strandedness check. 
    // We want to avoid adding another input stream here, so the path string of params.strandedness_gtf must be used.
    def gtf_path = gtf

    // Extract just the filename from the string to use for the local unzipped copy
    // We use .split('/')[-1] to get the filename from the full path string
    def gtf_filename = gtf.getName()
    
    // Determine the local unzipped name
    def gtf_local = gtf_filename.endsWith('.gz') ? gtf_filename.replace('.gz', '') : gtf_filename
    
    // Build the unzip command (streaming from the external path to the local work dir)
    def unzip_cmd = gtf_filename.endsWith('.gz') ? "gunzip -c ${gtf_path} > ${gtf_local}" : ""
    
    // If it wasn't gzipped, we use the original path directly; otherwise, use the local unzipped file
    def gtf_file = gtf_filename.endsWith('.gz') ? gtf_local : gtf_path


    """

    ${unzip_cmd}

    check_strandedness -g ${gtf_file} -k ${params.strandedness_ref} -r1 ${reads[0]} ${paired} > ${sampleID}_strandedness.txt 2>&1

    if grep -q "Data is likely" ${sampleID}_strandedness.txt; then
        
        data_type=`grep "Data is likely" ${sampleID}_strandedness.txt`
        
        if  [[ \$data_type == *RF* ]] ; then
            STRAND='reverse_stranded'
        elif [[ \$data_type == *FR* ]] ; then
            STRAND='forward_stranded'
        elif [[ \$data_type == *unstranded* ]] ; then
            STRAND='non_stranded'
        else
            echo "RNA Seq data does not fall into a likely stranded (max percent explained > 0.9) or unstranded layout (max percent explained < 0.6). Please check your data for low quality and contaminating reads before proceeding."; exit 42;
        fi
    
    else
    
        if [[ ${params.strandedness} == "reverse_stranded" || ${params.strandedness} == "forward_stranded" || ${params.strandedness} == "non_stranded" ]] ; then
            STRAND='${params.strandedness}'
        else
            echo "RNA Seq data does not fall into a likely stranded (max percent explained > 0.9) or unstranded layout (max percent explained < 0.6), and the parameter to override this: '--strandedness' is not set. Please check your data for low quality and contaminating reads"; exit 42;
        fi
    
    fi
    """
}

// Data is likely RF/fr-firststrand
// Data is likely FR/fr-secondstrand
// Data is likely unstranded
/*

# Script logic: 

if single_strand:
    fwd = float(result.iloc[2,0].replace('Fraction of reads explained by "++,--": ', ''))
    rev = float(result.iloc[3,0].replace('Fraction of reads explained by "+-,-+": ', ''))
else:
    fwd = float(result.iloc[2,0].replace('Fraction of reads explained by "1++,1--,2+-,2-+": ', ''))
    rev = float(result.iloc[3,0].replace('Fraction of reads explained by "1+-,1-+,2++,2--": ', ''))
if float(result.iloc[1,0].replace('Fraction of reads failed to determine: ', '')) > 0.50:
    print('Failed to determine strandedness of > 50% of reads.')
    print('If this is unexpected, try running again with a higher --nreads value')
if fwd_percent > 0.9:
    if single_strand:
        print('Over 90% of reads explained by "++,--"')
        print('Data is likely FR/fr-stranded')
    else:
        print('Over 90% of reads explained by "1++,1--,2+-,2-+"')
        print('Data is likely FR/fr-secondstrand')
elif rev_percent > 0.9:
    if single_strand:
        print('Over 90% of reads explained by "+-,-+"')
        print('Data is likely RF/rf-stranded')
    else:
        print('Over 90% of reads explained by "1+-,1-+,2++,2--"')
        print('Data is likely RF/fr-firststrand')
elif max(fwd_percent, rev_percent) < 0.6:
    print('Under 60% of reads explained by one direction')
    print('Data is likely unstranded')
else:
    print('Data does not fall into a likely stranded (max percent explained > 0.9) or unstranded layout (max percent explained < 0.6)')
    print('Please check your data for low quality and contaminating reads before proceeding')
*/
