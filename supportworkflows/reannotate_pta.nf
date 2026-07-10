#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {ANNOTATE_SV_WITH_CNV as ANNOTATE_SV_WITH_CNV_HUMAN;
        ANNOTATE_SV_WITH_CNV as ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL_HUMAN} from "../modules/r/annotate_sv_with_cnv"
include {FILTER_BEDPE as FILTER_BEDPE_HUMAN;
        FILTER_BEDPE as FILTER_BEDPE_SUPPLEMENTAL_HUMAN} from "../modules/r/filter_bedpe"

include {ANNOTATE_SV_WITH_CNV as ANNOTATE_SV_WITH_CNV_MOUSE;
        ANNOTATE_SV_WITH_CNV as ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL_MOUSE} from "../modules/r/annotate_sv_with_cnv_mouse"
include {FILTER_BEDPE as FILTER_BEDPE_MOUSE;
        FILTER_BEDPE as FILTER_BEDPE_SUPPLEMENTAL_MOUSE} from "../modules/r/filter_bedpe_mouse"

workflow REANNOTATE_PTA {
    if (!params.csv_input) {
        exit 1, "No input CSV file was specified with `--csv_input`. A CSV manifest is required. See the GitHub wiki (https://github.com/TheJacksonLaboratory/jds-nf-workflows/wiki/PTA-Pipeline-ReadMe) for information."
    }

    if (params.csv_input) {
        input_files = extract_csv(file(params.csv_input, checkIfExists: true))
        main_files = input_files.map{ it[0] } // main
        supp_files = input_files.map{ it[1] } // supplemental
    }

    if (params.gen_org == 'human') {
        ANNOTATE_SV_WITH_CNV_HUMAN(main_files, "main")
        ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL_HUMAN(supp_files, "supplemental")

        FILTER_BEDPE_HUMAN(ANNOTATE_SV_WITH_CNV_HUMAN.out.sv_genes_cnv_bedpe, "main")
        FILTER_BEDPE_SUPPLEMENTAL_HUMAN(ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL_HUMAN.out.sv_genes_cnv_bedpe, "supplemental")
    }
            
    if (params.gen_org == 'mouse') {
        ANNOTATE_SV_WITH_CNV_MOUSE(main_files, "main")
        ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL_MOUSE(supp_files, "supplemental")
        
        FILTER_BEDPE_MOUSE(ANNOTATE_SV_WITH_CNV_MOUSE.out.sv_genes_cnv_bedpe, "main")
        FILTER_BEDPE_SUPPLEMENTAL_MOUSE(ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL_MOUSE.out.sv_genes_cnv_bedpe, "supplemental")
   
    }
}


/*
These are the post annotation, and CNV annotated files. However, these must be reannotated with CNV.

# Human:
manta_gridss_sv_annotated_genes_cnv.bedpe
manta_gridss_sv_annotated_genes_cnv_supplemental.bedpe

# Mouse: 
manta_lumpy_delly_svaba_sv_annotated_genes_cnv.bedpe
manta_lumpy_delly_svaba_sv_annotated_genes_cnv_supplemental.bedpe

----

These are the CNV files that are post annotation. They are the same for mouse and human:
cnv_annotated_final.bed
cnv_annotated_supplemental.bed
*/


def extract_csv(csv_file) {
    def ANSI_RED = "\u001B[31m"
    def ANSI_RESET = "\u001B[0m"

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def sampleSheetLines = reader.readLines()
        def numberOfLinesInSampleSheet = sampleSheetLines.size()
        if (numberOfLinesInSampleSheet < 2) {
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, and at least one sample." + ANSI_RESET)
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }

        // Reopen the file to read and check the header
        def headerLine
        file(csv_file).withReader('UTF-8') { headerReader ->
            headerLine = headerReader.readLine()
        }
        def headers = headerLine.split(',').collect { it.trim() }
        def requiredHeaders = ['sampleID', 'sv_annotated_genes_cnv', 'annotated_genes_cnv_supplemental', 'cnv_annotated', 'cnv_annotated_supplemental']

        def requiredHeadersStr = requiredHeaders.collect { "'${it}'" }.join(', ')
  
        def missingHeaders = requiredHeaders.findAll { !headers.contains(it) }
        if (missingHeaders) {
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "Missing required header(s) in CSV file: ${missingHeaders.join(', ')}" + ANSI_RESET)
            System.err.println(ANSI_RED + "The csv file must have fields: ${requiredHeadersStr}" + ANSI_RESET)
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }
    }

    channel.from(csv_file).splitCsv(header: true)
        .map{ row ->

            // List of files to check for existence.
            def filesToCheck = [
                row.cnv_annotated,
                row.sv_annotated_genes_cnv,
                row.cnv_annotated_supplemental,
                row.annotated_genes_cnv_supplemental
            ]

            filesToCheck.each { filePath ->
                try {
                    file(filePath, checkIfExists: true)
                } catch (Exception e) {
                    System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                    def headerName = row.find { k, v -> v == filePath }?.key ?: "unknown"
                    System.err.println(ANSI_RED + "The file for header '${headerName}': " + filePath + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
                    System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                    System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                    System.exit(1)
                }
            }
            
            def main = [row.sampleID.toString(), 'normal_placeholder', 'tumor_placeholder', file(row.cnv_annotated, checkIfExists: true), file(row.sv_annotated_genes_cnv, checkIfExists: true)]
            def supp = [row.sampleID.toString(), 'normal_placeholder', 'tumor_placeholder', file(row.cnv_annotated_supplemental, checkIfExists: true), file(row.annotated_genes_cnv_supplemental, checkIfExists: true)]
            // Note: normal and tumor names are placeholders here, as they are not used in the reannotation steps. The module input tuples require them, but do not use the strings. 

            return [main, supp]

        } // end row map
    // end of channel.from, no bracket needed.
} // end extract_csv
