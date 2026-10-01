def extract_csv(csv_file) {
    file(csv_file).withReader('UTF-8') { reader ->
        def lines = reader.readLines()

        if (lines.size() < 2) {
            error "The CSV must contain a header and at least one sample."
        }

        def headers = lines.first().split(',').collect { header -> header.trim() }
        def requiredHeaders = ['SampleID', 'Message']
        def missingHeaders = requiredHeaders.findAll { header -> !headers.contains(header) }

        if (missingHeaders) {
            error "Missing required CSV header(s): ${missingHeaders.join(', ')}"
        }
    }

    channel.from(csv_file)
        .splitCsv(header: true)
        .map { row ->
            if (!row.SampleID || !row.Message) {
                error "SampleID and Message must not be empty: ${row}"
            }

            tuple(row.SampleID.toString(), row.Message.toString())
        }
}