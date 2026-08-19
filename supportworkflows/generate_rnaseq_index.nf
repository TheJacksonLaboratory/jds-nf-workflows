#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "../bin/help/generate_rnaseq_index.nf"
include {param_log} from "../bin/log/generate_rnaseq_index.nf"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {GUNZIP} from "../modules/utility_modules/gunzip"
include {AGAT_GFFTOGTF} from "../modules/agat/agat_gfftogtf"
include {GFFREAD_GFF3TOGTF} from "../modules/gffread/gffread_gff3togtf"
include {MODIFY_MGI_GTF} from "../modules/utility_modules/modify_mgi_gtf"
include {MAKE_CUSTOM_TRANSCRIPTOME} from "../modules/utility_modules/make_custom_transcriptome"
include {RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_BOWTIE2;
         RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_STAR} from "../modules/rsem/rsem_preparereference"
include {PICARD_CREATESEQUENCEDICTIONARY} from "../modules/picard/picard_createsequencedictionary"
include {UCSC_GTFTOGENEPRED} from "../modules/ucsc/ucsc_gtftogenepred"
include {GENERATE_RRNA_INTERVALS} from "../modules/utility_modules/generate_rrna_intervals"
include {KALLISTO_INDEX} from "../modules/kallisto/kallisto_index"


workflow GENERATE_RNASEQ_INDEX {

    // help if needed
    if (params.help){
        help()
        exit 0
    }

    // log params
    message = param_log()

    // Save params to a file for record-keeping
    workflow.onComplete {
        final_run_report(message)
    }

    if( !params.ref_gff3) {
    if (!params.ref_gtf && !params.ref_gff) {
        log.error "Either `--ref_gtf` or `--ref_gff` must be provided."
        exit 1
    }

    if ((params.ref_gtf != null && params.ref_gtf != '') && (params.ref_gff != null && params.ref_gff != '')) {
        log.error "Both `--ref_gtf` and `--ref_gff` were provided. Please provide only one."
        exit 1
    }
    if (params.ref_gff) {
        checkFileExists(params.ref_gff, "ref_gff")
    }

    } else {
        if (!params.ref_gtf && !params.ref_gff3) {
        log.error "Either `--ref_gtf` or `--ref_gff3` must be provided."
        exit 1
        }
        
        if ((params.ref_gtf != null && params.ref_gtf != '') && (params.ref_gff3 != null && params.ref_gff3 != '')) {
        log.error "Both `--ref_gtf` and `--ref_gff3` were provided. Please provide only one."
        exit 1
    }
    if (params.ref_gff3) {
        checkFileExists(params.ref_gff3, "ref_gff3")
    }

    }     

    if (params.ref_gtf) {
        checkFileExists(params.ref_gtf, "ref_gtf")
    }

    checkFileExists(params.ref_fa, "ref_fa")

    if (params.custom_gene_fasta) {
        checkFileExists(params.custom_gene_fasta, "custom_gene_fasta")
    }
    
    star_read_lengths = params.star_read_lengths instanceof String 
        ? channel.from(params.star_read_lengths.split(',').collect{ it -> it.trim().toInteger()})
        : channel.from(params.star_read_lengths.collect{ it -> it.toInteger()})

    if (params.ref_gff) {
        AGAT_GFFTOGTF(params.ref_gff)
        proc_gtf = AGAT_GFFTOGTF.out.gtf
    } else if (params.ref_gff3) {
        if (params.ref_gff3.endsWith('.gz')) {
            GUNZIP(params.ref_gff3)
            GFFREAD_GFF3TOGTF(GUNZIP.out.gunzip_file)
            proc_gtf = GFFREAD_GFF3TOGTF.out.gtf
        } else {
            GFFREAD_GFF3TOGTF(params.ref_gff3)
            proc_gtf = GFFREAD_GFF3TOGTF.out.gtf
        }
    } else {
        if (params.ref_gtf.endsWith('.gz')) {
            GUNZIP(params.ref_gtf)
            proc_gtf = GUNZIP.out.gunzip_file
        } else {
            proc_gtf = channel.fromPath(params.ref_gtf)
        }
    }

    if (params.mgi) {
        MODIFY_MGI_GTF(proc_gtf, params.ref_table)
        proc_gtf = MODIFY_MGI_GTF.out.gtf
    }

    if (params.custom_gene_fasta) {
    
        make_custom_input = proc_gtf
                            .map{it -> [params.ref_fa, it, params.custom_gene_fasta]}

        MAKE_CUSTOM_TRANSCRIPTOME(make_custom_input)

        bowtie2_input = MAKE_CUSTOM_TRANSCRIPTOME.out.concat_fasta.combine(MAKE_CUSTOM_TRANSCRIPTOME.out.concat_gtf).map{it -> [it[0], it[1], 'bowtie2', '']}
        star_build_set = MAKE_CUSTOM_TRANSCRIPTOME.out.concat_fasta.combine(MAKE_CUSTOM_TRANSCRIPTOME.out.concat_gtf).combine(star_read_lengths).map{it -> [it[0], it[1], 'STAR', it[2]]}
        
        fasta = MAKE_CUSTOM_TRANSCRIPTOME.out.concat_fasta
        gtf = MAKE_CUSTOM_TRANSCRIPTOME.out.concat_gtf
    
    } else {
        bowtie2_input = proc_gtf
        .map{ it -> [params.ref_fa, it, 'bowtie2', ''] }

        star_build_set = proc_gtf
        .map{ it -> [params.ref_fa, it, 'star'] }.combine(star_read_lengths)
        
        fasta = params.ref_fa
        gtf = proc_gtf
    }

    RSEM_PREPAREREFERENCE_BOWTIE2(bowtie2_input)
    
    RSEM_PREPAREREFERENCE_STAR(star_build_set)

    PICARD_CREATESEQUENCEDICTIONARY(fasta)

    GENERATE_RRNA_INTERVALS(gtf, PICARD_CREATESEQUENCEDICTIONARY.out.dict)
    
    UCSC_GTFTOGENEPRED(gtf)

    KALLISTO_INDEX(RSEM_PREPAREREFERENCE_BOWTIE2.out.transcripts)

}

def checkFileExists(filePath, name) {
    if (filePath && !file(filePath).exists()) {
        log.error "File not found: ${filePath} (${name})"
        exit 1
    }
}

// Workflow adapted from: https://github.com/KU-GDSC/workflows
