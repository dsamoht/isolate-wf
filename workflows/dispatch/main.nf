#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow DISPATCH {

    main:
    def samplesheet_file = file(params.input, checkIfExists: true)
    def lines = samplesheet_file.text.readLines()
    lines.eachWithIndex { line, idx ->
        def columns = line.split(',', -1)
        if (columns.size() != 4) {
            def line_num = idx
            exit 1, "sample sheet error : ${columns.size()} field(s) found at line ${line_num}; exactly 4 fields (3 commas) are expected."
        }
    }

    ch_samplesheet = channel
        .fromPath(samplesheet_file)
        .ifEmpty { exit 1, "sample sheet not found: `${params.input}`" }
        .splitCsv(header: true)
        .map { row ->
            def required = [
                'sample_id',
                'long_reads',
                'short_reads_1',
                'short_reads_2'
            ]

            required.each { col ->
                if (!row.containsKey(col)) {
                    exit 1, "sample sheet is missing required column: `${col}`"
                }
            }

            if (!row.sample_id) {
                exit 1, "sample sheet row has empty `sample_id`"
            }

            if (!row.long_reads) {
                exit 1, "sample sheet row for `${row.sample_id}` has empty `long_reads` (required for `autocycler`)"
            }
            
            def has_sr1 = row.short_reads_1 ? true : false
            def has_sr2 = row.short_reads_2 ? true : false
            if (has_sr1 != has_sr2) {
                exit 1, "sample `${row.sample_id}` has unpaired short reads - both `short_reads_1` and `short_reads_2` are required"
            }

            def long_reads   = file(row.long_reads, checkIfExists: true)
            def short_reads1 = row.short_reads_1  ? file(row.short_reads_1, checkIfExists: true)  : false
            def short_reads2 = row.short_reads_2  ? file(row.short_reads_2, checkIfExists: true)  : false

            return [
                sample_id     : row.sample_id,
                long_reads    : long_reads,
                sr1           : short_reads1,
                sr2           : short_reads2
            ]
        }

    // Validate sample sheet to ensure no duplicate sample IDs
    ch_validate = ch_samplesheet
        .map { row -> [ row.sample_id, row ] }
        .groupTuple(by: 0)
        .map { sample_id, rows ->
            if (rows.size() > 1) {
                exit 1, "sample sheet error: duplicate sample_id found: '${sample_id}'"
            }
        }
        .collect()

    emit:
    samplesheet = ch_samplesheet
    validate    = ch_validate

}
