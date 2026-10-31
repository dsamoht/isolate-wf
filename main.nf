#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dsamoht/isolate-wf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/dsamoht/isolate-wf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AUTOCYCLER_ASSEMBLE         } from './modules/local/autocycler/assemble'
include { AUTOCYCLER_CLUSTER          } from './modules/local/autocycler/cluster'
include { AUTOCYCLER_COMBINE          } from './modules/local/autocycler/combine'
include { AUTOCYCLER_COMPRESS         } from './modules/local/autocycler/compress'
include { AUTOCYCLER_GENOME_SIZE      } from './modules/local/autocycler/genome_size'
include { AUTOCYCLER_SUBSAMPLE        } from './modules/local/autocycler/subsample'
include { AUTOCYCLER_TRIM_AND_RESOLVE } from './modules/local/autocycler/trim_and_resolve'
include { CHECKM2                     } from './modules/local/checkm2'
include { DNAAPLER                    } from './modules/local/dnaapler'
include { GENOMAD                     } from './modules/local/genomad'
include { LONGREAD_QC                 } from './workflows/longread_qc'
include { MULTIQC                     } from './modules/local/multiqc'
include { PIPELINE_INITIALISATION     } from './subworkflows/local/pipeline_initialisation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BOOLEAN PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Declared with their type so that Nextflow converts `--skip_qc false` on the command line
// to false: undeclared, the value stays the string "false", which Groovy reads as true.
// Defaults are set in nextflow.config. `--help` is left out, it also takes a string.
params {
    skip_qc          : Boolean
    skip_nanoplot    : Boolean
    skip_porechop    : Boolean
    skip_flye        : Boolean
    skip_plassembler : Boolean
    skip_canu        : Boolean
    skip_metamdbg    : Boolean
    skip_miniasm     : Boolean
    skip_myloasm     : Boolean
    skip_raven       : Boolean
    skip_necat       : Boolean
    skip_nextdenovo  : Boolean
    skip_dnaapler    : Boolean
    skip_checkm2     : Boolean
    skip_genomad     : Boolean
    skip_multiqc     : Boolean
    monochrome_logs  : Boolean
    help_full        : Boolean
    show_hidden      : Boolean
    version          : Boolean
    validate_params  : Boolean
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISOLATE_WF {

    take:
    ch_samplesheet // channel: sample rows read from the sample sheet

    main:
    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    // Channel of input long reads
    ch_raw_long_reads = ch_samplesheet
        .map { sample -> [ [ id: sample.sample_id ], sample.long_reads ] }

    //
    // SUBWORKFLOW: long read quality control
    //
    ch_qc_long_reads_out = LONGREAD_QC(ch_raw_long_reads)
    ch_qc_long_reads     = ch_qc_long_reads_out.long_reads_qc
    ch_versions          = ch_versions.mix(ch_qc_long_reads_out.versions)
    ch_multiqc_files     = ch_multiqc_files.mix(ch_qc_long_reads_out.multiqc_files)

    //
    // MODULE: estimate the genome size of each sample
    //
    ch_genome_size = AUTOCYCLER_GENOME_SIZE(ch_qc_long_reads).size
    ch_versions = ch_versions.mix(AUTOCYCLER_GENOME_SIZE.out.versions.first())

    //
    // MODULE: split the reads into subsets, one per assembly
    //
    ch_subsampled_reads = AUTOCYCLER_SUBSAMPLE(ch_qc_long_reads.join(ch_genome_size)).reads
    ch_subsample_yaml   = AUTOCYCLER_SUBSAMPLE.out.yaml
    ch_versions = ch_versions.mix(AUTOCYCLER_SUBSAMPLE.out.versions.first())

    // Every read subset is assembled with every assembler not turned off with `--skip_<assembler>`
    def assemblers = [
        'flye', 'plassembler', 'canu', 'metamdbg',
        'miniasm', 'myloasm', 'raven', 'necat', 'nextdenovo'
    ].findAll { assembler -> !params["skip_${assembler}"] }

    // Flatten subsamples so each read subset is assembled independently. Each task carries
    // the number of assemblies its sample expects, so a sample is regrouped as soon as its
    // own assemblies are done instead of waiting for every sample
    ch_assembly_inputs = ch_subsampled_reads
        .join(ch_genome_size)
        .flatMap { meta, read_files, size ->
            def subsamples = [read_files].flatten()
            def n_assemblies = subsamples.size() * assemblers.size()
            subsamples.collectMany { read_file ->
                def subsample_id = read_file.baseName.tokenize('_').last()
                assemblers.collect { assembler ->
                    def asm_meta = [
                        id           : "${meta.id}_subsample_${subsample_id}_${assembler}",
                        sample       : meta.id,
                        assembler    : assembler,
                        n_assemblies : n_assemblies
                    ]
                    [ asm_meta, read_file, size, params.read_type ]
                }
            }
        }

    //
    // MODULE: assemble each read subset with each assembler, applying the clustering
    // and consensus weights
    //
    ch_raw_assemblies = AUTOCYCLER_ASSEMBLE(ch_assembly_inputs).fasta
    ch_versions = ch_versions.mix(AUTOCYCLER_ASSEMBLE.out.versions.first())

    // Regroup assemblies by sample; a failed assembler contributes an empty list
    ch_all_assemblies = ch_raw_assemblies
        .map { meta, fastas -> [ groupKey(meta.sample, meta.n_assemblies), fastas ] }
        .groupTuple()
        .map { sample, fastas -> [ [ id: sample.toString() ], fastas.flatten() ] }
        .branch { _meta, fastas ->
            assembled: fastas
            failed:    true
        }

    ch_all_assemblies.failed.subscribe { meta, _fastas ->
        log.warn("Sample '${meta.id}': every assembler failed, no final assembly")
    }

    //
    // MODULE: compress the input assemblies into a unitig graph
    //
    ch_compressed = AUTOCYCLER_COMPRESS(ch_all_assemblies.assembled).graph
    ch_versions = ch_versions.mix(AUTOCYCLER_COMPRESS.out.versions.first())

    //
    // MODULE: cluster the input contigs into putative genomic sequences
    //
    ch_clustered = AUTOCYCLER_CLUSTER(ch_compressed)
    ch_versions = ch_versions.mix(AUTOCYCLER_CLUSTER.out.versions.first())

    ch_clustered.autocycler_out
        .join(ch_clustered.qc_clusters, remainder: true)
        .filter { _meta, _autocycler_out, clusters -> clusters == null }
        .subscribe { meta, _autocycler_out, _clusters ->
            log.warn("Sample '${meta.id}': no cluster passed Autocycler QC, no final assembly")
        }

    // Extract individual QC-pass clusters (e.g., cluster_1, cluster_2) for parallel
    // processing, keeping the cluster count for the regroup after trimming and resolving
    ch_qc_clusters = ch_clustered.qc_clusters
        .flatMap { meta, cluster_dirs ->
            def clusters = [cluster_dirs].flatten()
            clusters.collect { cluster_dir ->
                [ [ id: "${meta.id}_${cluster_dir.baseName}", sample: meta.id, n_clusters: clusters.size() ], cluster_dir ]
            }
        }

    //
    // MODULE: trim and resolve each QC-pass cluster
    //
    ch_resolved_clusters = AUTOCYCLER_TRIM_AND_RESOLVE(ch_qc_clusters).cluster
    ch_versions = ch_versions.mix(AUTOCYCLER_TRIM_AND_RESOLVE.out.versions.first())

    // Collect all resolved cluster directories to combine per sample
    ch_all_resolved = ch_resolved_clusters
        .map { meta, cluster_dir -> [ groupKey(meta.sample, meta.n_clusters), cluster_dir ] }
        .groupTuple()
        .map { sample, cluster_dirs -> [ [ id: sample.toString() ], cluster_dirs ] }

    //
    // MODULE: combine the resolved clusters into a final assembly and tabulate the metrics
    //
    AUTOCYCLER_COMBINE(
        ch_all_resolved
            .join(ch_clustered.autocycler_out)
            .join(ch_subsample_yaml)
    )
    ch_final_assembly      = AUTOCYCLER_COMBINE.out.fasta
    ch_consensus_graph     = AUTOCYCLER_COMBINE.out.gfa
    ch_autocycler_metrics  = AUTOCYCLER_COMBINE.out.metrics
    ch_versions = ch_versions.mix(AUTOCYCLER_COMBINE.out.versions.first())

    // One archive of the Autocycler directory per sample. The one from AUTOCYCLER_COMBINE
    // also holds the trimmed and resolved clusters and the consensus; a sample with no
    // cluster passing QC never gets there and falls back on the one from AUTOCYCLER_CLUSTER
    ch_autocycler_archive = ch_clustered.archive
        .join(AUTOCYCLER_COMBINE.out.archive, remainder: true)
        .map { meta, cluster_archive, combine_archive -> [ meta, combine_archive ?: cluster_archive ] }

    // One table for the whole run: the header of the first file, then a line per sample.
    // Sorted on the file name, which starts with the sample id, so the order depends neither
    // on which sample finished first nor on the work directory hashes
    ch_autocycler_metrics_summary = ch_autocycler_metrics
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'autocycler_metrics.tsv', keepHeader: true, sort: { tsv -> tsv.name })

    //
    // MODULE: reorient circular contigs to start at dnaA (or repA/terL)
    //
    ch_reoriented_assembly = channel.empty()
    if (!params.skip_dnaapler) {
        ch_reoriented_assembly = DNAAPLER(ch_final_assembly).fasta
        ch_versions = ch_versions.mix(DNAAPLER.out.versions.first())
    }

    //
    // MODULE: estimate completeness and contamination of the final assembly
    //
    ch_checkm2_report = channel.empty()
    if (!params.skip_checkm2) {
        ch_checkm2_report = CHECKM2(ch_final_assembly).report
        ch_versions       = ch_versions.mix(CHECKM2.out.versions.first())
        ch_multiqc_files  = ch_multiqc_files.mix(ch_checkm2_report.map { _meta, report -> report })
    }

    //
    // MODULE: identify plasmids, viruses and proviruses in the final assembly
    //
    // geNomad has no default database to fall back on: without `--genomad_db` it is skipped
    // as if `--skip_genomad` had been given
    ch_genomad_results = channel.empty()
    if (!params.skip_genomad) {
        if (params.genomad_db) {
            ch_genomad_results = GENOMAD(ch_final_assembly, file(params.genomad_db, checkIfExists: true)).results
            ch_versions        = ch_versions.mix(GENOMAD.out.versions.first())
        } else {
            log.info("No --genomad_db given, skipping geNomad")
        }
    }

    //
    // Collate software versions
    //
    // Collect the contents rather than the paths: `sort: true` orders the entries by their
    // value, and a channel of paths sorts by work directory, which is a fresh random hash on
    // every run. Reading each file first makes both the deduplication and the order depend on
    // the versions themselves, so the published file is the same from one run to the next.
    ch_collated_versions = ch_versions
        .map { versions -> versions.text }
        .unique()
        .collectFile(name: 'isolate-wf_software_mqc_versions.yml', sort: true)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_config = channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config
            ? channel.fromPath(params.multiqc_config, checkIfExists: true)
            : channel.empty()
        ch_multiqc_logo = params.multiqc_logo
            ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
            : channel.empty()

        // MultiQC writes no report when its search turns up nothing a module can parse. The
        // software versions file alone does not count: it only fills the versions table and
        // renders no section. `collect()` emitting nothing therefore means there is no report
        // to write — with `--skip_qc --skip_checkm2`, for instance — and MULTIQC is not
        // submitted at all instead of failing on a report it never wrote.
        ch_multiqc_input = ch_multiqc_files
            .collect()
            .combine(ch_collated_versions)

        MULTIQC(
            ch_multiqc_input,
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    versions            = ch_collated_versions
    multiqc_report      = ch_multiqc_report
    nanoplot_raw_html   = ch_qc_long_reads_out.nanoplot_raw_html
    nanoplot_qc_html    = ch_qc_long_reads_out.nanoplot_qc_html
    porechop_log        = ch_qc_long_reads_out.porechop_log
    qc_long_reads       = ch_qc_long_reads
    final_assembly             = ch_final_assembly
    consensus_graph            = ch_consensus_graph
    autocycler_metrics         = ch_autocycler_metrics
    autocycler_metrics_summary = ch_autocycler_metrics_summary
    clustering                 = ch_clustered.clustering
    autocycler_archive         = ch_autocycler_archive
    reoriented_assembly        = ch_reoriented_assembly
    checkm2_report             = ch_checkm2_report
    genomad_results            = ch_genomad_results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN THE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: validate the parameters and read the sample sheet
    //
    PIPELINE_INITIALISATION(
        params.validate_params,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        params.version
    )

    //
    // WORKFLOW: run the main analysis
    //
    ISOLATE_WF(PIPELINE_INITIALISATION.out.samplesheet)

    publish:
    versions            = ISOLATE_WF.out.versions
    multiqc_report      = ISOLATE_WF.out.multiqc_report
    nanoplot_raw_html   = ISOLATE_WF.out.nanoplot_raw_html
    nanoplot_qc_html    = ISOLATE_WF.out.nanoplot_qc_html
    porechop_log        = ISOLATE_WF.out.porechop_log
    qc_long_reads       = ISOLATE_WF.out.qc_long_reads
    final_assembly             = ISOLATE_WF.out.final_assembly
    consensus_graph            = ISOLATE_WF.out.consensus_graph
    autocycler_metrics         = ISOLATE_WF.out.autocycler_metrics
    autocycler_metrics_summary = ISOLATE_WF.out.autocycler_metrics_summary
    clustering                 = ISOLATE_WF.out.clustering
    autocycler_archive         = ISOLATE_WF.out.autocycler_archive
    reoriented_assembly        = ISOLATE_WF.out.reoriented_assembly
    checkm2_report             = ISOLATE_WF.out.checkm2_report
    genomad_results            = ISOLATE_WF.out.genomad_results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PUBLISH TARGETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    qc_long_reads {
        path { meta, _file -> "sample_${meta.id}/reads/post_qc" }
    }
    nanoplot_raw_html {
        path { meta, _file -> "sample_${meta.id}/quality_assessment/nanoplot/raw" }
    }
    nanoplot_qc_html {
        path { meta, _file -> "sample_${meta.id}/quality_assessment/nanoplot/post_qc" }
    }
    porechop_log {
        path { meta, _file -> "sample_${meta.id}/quality_control/porechop" }
    }
    final_assembly {
        path { meta, _file -> "sample_${meta.id}/autocycler_final" }
    }
    consensus_graph {
        path { meta, _file -> "sample_${meta.id}/autocycler_final" }
    }
    autocycler_metrics {
        path { meta, _file -> "sample_${meta.id}/autocycler_final" }
    }
    autocycler_metrics_summary {
        path "."
    }
    clustering {
        path { meta, _files -> "sample_${meta.id}/autocycler_clustering" }
    }
    autocycler_archive {
        path { meta, _file -> "sample_${meta.id}/autocycler_clustering" }
    }
    reoriented_assembly {
        path { meta, _file -> "sample_${meta.id}/dnaapler" }
    }
    checkm2_report {
        path { meta, _file -> "sample_${meta.id}/checkm2" }
    }
    genomad_results {
        path { meta, dir -> dir >> "sample_${meta.id}/genomad" }
    }
    multiqc_report {
        path "multiqc"
    }
    versions {
        path "pipeline_info"
    }
}
