//
// Subworkflow with functionality specific to the dsamoht/isolate-wf pipeline
//

include { paramsHelp         } from 'plugin/nf-schema'
include { paramsSummaryLog   } from 'plugin/nf-schema'
include { samplesheetToList  } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    validate_params // boolean: validate parameters against nextflow_schema.json
    input           // string : path to the sample sheet
    help            // boolean or string: display the help message
    help_full       // boolean: display the full help message
    show_hidden     // boolean: display hidden parameters in the help message
    version         // boolean: display the pipeline version

    main:

    //
    // Banner: same ASCII logo and description as the README, printed on top of the help
    // message and of every run
    //
    log.info(logo())

    //
    // Print the help message built from nextflow_schema.json and exit
    //
    if (help || help_full) {
        log.info(paramsHelp(
            [
                command   : "nextflow run dsamoht/isolate-wf -profile <docker/singularity/apptainer> --input samplesheet.csv --outdir <OUTDIR>",
                showHidden: show_hidden,
                fullHelp  : help_full,
            ],
            help instanceof String && help != "true" ? help : "",
        ))
        removeEmptyWorkDir()
        exit(0)
    }

    //
    // Print the pipeline version and exit
    //
    if (version) {
        log.info("${workflow.manifest.name} v${workflow.manifest.version}")
        removeEmptyWorkDir()
        exit(0)
    }

    //
    // Validate parameters and print a summary of the ones that differ from the defaults
    //
    if (validate_params) {
        validateParameters()
    }

    // `outputDir` follows `params.outdir` (nextflow.config), and Nextflow falls back to
    // publishing in `./results` when it is null. Parameter validation already requires
    // `--outdir`, but it can be turned off, so the run is stopped here as well: the results
    // are always written to the directory asked for on the command line, never anywhere else.
    if (!params.outdir) {
        error("No output directory: pass --outdir. Without it the results would be published to './results' instead.")
    }

    // Autocycler builds its consensus from the assemblies of the assemblers left on, so
    // turning every one of them off would leave it nothing to work with.
    def assemblers = ['flye', 'plassembler', 'canu', 'metamdbg', 'miniasm', 'myloasm', 'raven', 'necat', 'nextdenovo']
    if (assemblers.every { assembler -> params["skip_${assembler}"] }) {
        error("Every assembler is skipped: leave at least one of ${assemblers.collect { assembler -> "--skip_${assembler}" }.join(', ')} off.")
    }

    log.info(paramsSummaryLog(workflow))

    //
    // Read the sample sheet. Per-row validation lives in assets/schema_input.json;
    // nf-schema returns [ meta, long_reads, short_reads_1, short_reads_2 ] with an empty
    // list for any column left blank.
    //
    def samples = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        .collect { meta, long_reads, short_reads_1, short_reads_2 ->
            [
                sample_id  : meta.id,
                long_reads : long_reads,
                sr1        : short_reads_1 ?: false,
                sr2        : short_reads_2 ?: false
            ]
        }

    //
    // Checks that span several rows of the sample sheet, which a per-row JSON schema
    // cannot express. They run on the parsed sample sheet, before the channel the rest of
    // the pipeline reads from exists, so a failure stops the run before any task is submitted.
    //
    validateSamples(samples)

    ch_samplesheet = channel.fromList(samples)

    emit:
    samplesheet = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// ASCII logo, pipeline description and repository details, taken from the manifest so
// that the banner cannot drift from the version being run.
//
def logo() {
    return """

 _           _       _                      __
(_)___  ___ | | __ _| |_ ___     __      __/ _|
| / __|/ _ \\| |/ _` | __/ _ \\____\\ \\ /\\ / / |_
| \\__ \\ (_) | | (_| | ||  __/_____\\ V  V /|  _|
|_|___/\\___/|_|\\__,_|\\__\\___|      \\_/\\_/ |_|

${workflow.manifest.description}

     Github: ${workflow.manifest.homePage}
     Version: v${workflow.manifest.version}

"""
}

//
// Nextflow creates the work directory when the session starts, before this script runs, so
// a run that only prints the help message or the version leaves an empty `work/` behind in
// whatever directory it was called from. Remove it again, but only if it is still empty,
// i.e. if this run is the one that created it.
//
def removeEmptyWorkDir() {
    if (workflow.workDir.exists() && workflow.workDir.list().size() == 0) {
        workflow.workDir.deleteDir()
    }
}

//
// Validate the sample sheet as a whole: sample ids name every output directory and file
// of a sample, so two rows sharing one would overwrite each other's results.
//
def validateSamples(rows) {
    def duplicate_ids = rows
        .countBy { row -> row.sample_id }
        .findAll { _id, count -> count > 1 }
        .keySet()
    if (duplicate_ids) {
        error("sample sheet has repeated sample ids: ${duplicate_ids.sort().join(', ')}. Sample ids must be unique.")
    }

    // The same read file must not be claimed by two samples
    def path_to_samples = rows.groupBy { row -> row.long_reads.toString() }
    path_to_samples.each { path, path_rows ->
        if (path_rows.size() > 1) {
            error("read file '${path}' is shared by multiple samples: ${path_rows.collect { row -> row.sample_id }.join(', ')}")
        }
    }
}
