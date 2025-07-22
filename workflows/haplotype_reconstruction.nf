#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Nextflow pipeline for sample QC and haplotype reconstruction
// on genetically diverse mice

// import modules
include {help} from "${projectDir}/bin/help/haplotype_reconstruction.nf"
include {param_log} from "${projectDir}/bin/log/haplotype_reconstruction.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv_qtl.nf"
include {GS_TO_QTL2} from "${projectDir}/modules/qtl2/geneseek2qtl2"
include {WRITE_CROSS} from "${projectDir}/modules/qtl2/write_cross"
include {GENOPROBS} from "${projectDir}/modules/qtl2/genoprobs"
include {CONCAT_GENOPROBS} from "${projectDir}/modules/qtl2/concat_genoprobs"
include {CONCAT_INTENSITIES} from "${projectDir}/modules/qtl2/concat_intensities"
include {UPDATE_FILES} from "${projectDir}/modules/qtl2/update_files"
include {QC_REPORT} from "${projectDir}/modules/r/render_QC_markdown"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// Make channel of consensus files (GigaMUGA)
founder_genos   = Channel.fromPath(params.gm_cc_do_founder_genotypes).collect()
gmaps           = Channel.fromPath(params.gm_gmaps).collect()
pmaps           = Channel.fromPath(params.gm_pmaps).collect()
consensus_files = founder_genos.concat(gmaps)
                               .concat(pmaps)
                               .flatten().collect()


// QC and Haplotype Reconstruction Workflow
workflow HAPLOTYPE_RECONSTRUCTION {

    project_ch = extract_csv(params.csv)

    if(params.rerun){
        // If rerun is true, we need to reprocess final report files
        // and re-run the pipeline
        println("FinalReport files supplied; performing haplotype reconstruction.")
        project_ch.map{it -> [  it[0], it[1].finalreport_file,
                                it[1].covar_file, it[1].cross_type]}.set{rerun_ch}
        
        // Process FinalReport File
        GS_TO_QTL2(rerun_ch)
        
        // Gather intensities by project id
        sexchr_intensities = GS_TO_QTL2.out.qtl2ints
                                    .groupTuple(by: 0)
                                    .map {it -> [it[0], it[1].flatten()]}

        all_intensities = GS_TO_QTL2.out.qtl2intsfst.groupTuple(by: 0)
        qc_intensities = sexchr_intensities.join(all_intensities, by: [0, 0])

        // Concatenate intensities across projects for QC
        CONCAT_INTENSITIES(qc_intensities)
        metadata = GS_TO_QTL2.out.qtl2meta.combine(CONCAT_INTENSITIES.out.dedup_samples, by: 0)
        sampleGenos = GS_TO_QTL2.out.sampleGenos

        // Write control file
        WRITE_CROSS(metadata, sampleGenos, consensus_files)

        // Initial haplotype reconstruction
        GENOPROBS(WRITE_CROSS.out.cross)

        // Gather genoprobs by project id
        project_genoprobs = GENOPROBS.out.genoprobs.groupTuple(by: 0)

        // Concatenate genoprobs across projects and perform marker QC
        CONCAT_GENOPROBS(project_genoprobs)

        // Join all the hr elements
        qc_data = CONCAT_GENOPROBS.out.qtl2_files.join(CONCAT_GENOPROBS.out.qc_files).join(CONCAT_INTENSITIES.out.concat_intensities)

        // Render the QC report
        QC_REPORT(qc_data)

    } else {
        // If rerun is false, we can jump to correcting ids or removing markers
        project_ch.map{it -> [  it[0], it[1].covar_file, it[1].cross_file,
                                it[1].genoprobs_file, it[1].alleleprobs_file,
                                it[1].viterbi_file, it[1].kinship_file, it[1].marker_file, it[1].cross_type ]}.set{update_ch}
        println("Haplotype reconstuctions supplied; updating files with supplied metadata.")
        println("Update individual ids: ${params.correct_ids}")
        println("Remove bad markers: ${params.remove_markers}")
        
        // Take probs, cross file, updated covar, marker files and update everything based on previous QC
        UPDATE_FILES(update_ch)
    }
}


    
