// Function to extract information (meta data + file(s)) from csv file(s)
// https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L1084

ANSI_RED = "\u001B[31m";
ANSI_RESET = "\u001B[0m";

def extract_csv(csv_file) {
    return Channel.fromPath(csv_file)
        .splitCsv(header: true)
        .map{ row ->
            if (!(row.project_id)){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "Missing project_id field in csv file header" + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }
            [row.project_id.toString(), row]
        }
        .map{ meta, row ->
            def meta_map = [:]
            meta_map.id = meta
            meta_map.finalreport_file = row.finalreport_file ?: 'NA'
            meta_map.covar_file = row.covar_file ?: 'NA'
            meta_map.cross_type = row.cross_type ?: 'NA'
            meta_map.marker_file = row.marker_file ?: 'NA'
            meta_map.genoprobs_file = row.genoprobs_file ?: 'NA'
            meta_map.alleleprobs_file = row.alleleprobs_file ?: 'NA'
            meta_map.cross_file = row.cross_file ?: 'NA'
            meta_map.viterbi_file = row.viterbi_file ?: 'NA'
            meta_map.kinship_file = row.kinship_file ?: 'NA'
            meta_map.pheno_file = row.pheno_file ?: 'NA'
            
            if (params.workflow == "sample_qc_haplotype_reconstructions" && meta_map.finalreport_file == 'NA') {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "`workflow = sample_qc_haplotype_reconstructions` specified but `finalreport_file` field is missing in the CSV manifest. Please add the `finalreport_file` field to the manifest and restart the run." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
            }

            if (params.workflow == "qtl_mapping" && meta_map.pheno_file == 'NA') {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "`workflow = qtl_mapping` specified but `pheno_file` field is missing in the CSV manifest. Please add the `pheno_file` field to the manifest and restart the run." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
            }

            // Check if `rerun` is specified and `finalreport_file` is missing
            if (params.rerun && meta_map.finalreport_file == 'NA') {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "`rerun` specified, but `finalreport_file` field is missing in the CSV manifest. Please add the `finalreport_file` field to the manifest and restart the run." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
            }

            // Add validation check for rerun parameter
            if (params.rerun && (meta_map.cross_type == 'NA' || meta_map.covar_file == 'NA')) {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "`rerun` specified, but `cross_type` or `covar_file` fields are missing in the CSV manifest. Please add these fields to the manifest and restart the run." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
            }
        
            // Add validation check for marker file parameter if rerun is NOT specified and remove_markers IS specified
            if (!(params.rerun) && params.remove_markers && meta_map.marker_file == 'NA') {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "`rerun` is false and `remove_markers` is true, but `marker_file` field is missing in the CSV manifest. Please add this field to the manifest and restart the run." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
            }

            if (!(params.rerun) && params.correct_ids && (meta_map.genoprobs_file == 'NA' || meta_map.alleleprobs_file == 'NA' || meta_map.cross_file == 'NA')) {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "`rerun` is false and `correct_ids` is true, but one or more required fields are missing in the CSV manifest." + ANSI_RESET)
            System.err.println(ANSI_RED + "Required fields: genoprobs_file, alleleprobs_file, cross_file" + ANSI_RESET)
            System.err.println(ANSI_RED + "Please add these fields to the manifest and restart the run." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
            }

            [meta_map.id, meta_map]
        }
}
