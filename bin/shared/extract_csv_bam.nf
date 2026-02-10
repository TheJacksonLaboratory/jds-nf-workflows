// Function to extract information (meta data + file(s)) from csv file(s)
// https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L1084

ANSI_RED = "\u001B[31m";
ANSI_RESET = "\u001B[0m";

def extract_csv_bam(csv_file) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
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
        def requiredHeaders = ['sampleID', 'bam', 'bai']

        if (params.deepvariant) {
            requiredHeaders << 'sex'
        }

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

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!(row.sampleID) | !(row.bam) | !(row.bai)){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "Missing row data in field: 'sampleID' or 'bam' or 'bai'. These fields can not be empty." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }
            [row.sampleID.toString(), row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes ->

        def meta = [:]

        // Meta data to identify samplesheet
        if (row.sampleID) meta.sampleID = row.sampleID.toString()

        // Parse optional "sex" field
        if (row.sex) meta.sex = row.sex.toString()
        else meta.sex = 'NA'
        
        // Define the ID field for the sample.
        meta.id = row.sampleID.toString()
        
        // defines the number of files for each sample. 
        meta.size = size

        try {
            file(row.bam, checkIfExists: true)
        }
        catch (Exception e) {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "The file: " + row.bam + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }
        try {
            file(row.bai, checkIfExists: true)
        }
        catch (Exception e) {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "The file: " + row.bai + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }


        return [meta.id, meta, row.bam, row.bai]

    }
}



def extract_csv_bam_rnaseq(csv_file) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
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
        def requiredHeaders = ['sampleID', 'bam']

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

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!(row.sampleID) | !(row.bam)){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "Missing row data in field: 'sampleID' or 'bam'. These fields can not be empty." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }
            [row.sampleID.toString(), row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes ->

        def meta = [:]

        // Meta data to identify samplesheet
        if (row.sampleID) meta.sampleID = row.sampleID.toString()

        // Define the ID field for the sample.
        meta.id = row.sampleID.toString()
        
        // defines the number of files for each sample. 
        meta.size = size

        try {
            file(row.bam, checkIfExists: true)
        }
        catch (Exception e) {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "The file: " + row.bam + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }

        return [meta.id, row.bam]

    }
}
