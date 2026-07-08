#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DISPATCH         } from './workflows/dispatch'
include { LONGREAD_QC      } from './workflows/longread_qc'
include { GENOME_SIZE      } from './modules/autocycler/genome_size'
include { SUBSAMPLE        } from './modules/autocycler/subsample'
include { ASSEMBLE         } from './modules/autocycler/assemble'
include { APPLY_WEIGHTS    } from './modules/autocycler/apply_weights'
include { COMPRESS         } from './modules/autocycler/compress'
include { CLUSTER          } from './modules/autocycler/cluster'
include { TRIM_AND_RESOLVE } from './modules/autocycler/trim_and_resolve'
include { COMBINE          } from './modules/autocycler/combine'

def helpMessage() {
    log.info """
    
nextflow wrapper for `autocycler` genome assembly pipeline

Usage:
     nextflow run main.nf -profile docker --input samplesheet.csv --outdir PATH
Input:
     --input         : path to input sample sheet (CSV format)
     --outdir        : path to output directory
     --read_type     : read type (ont_r9, ont_r10, pacbio_clr, pacbio_hifi) [default: ont_r10]
     --skip_qc       : skip long read qc/qa steps
     --skip_nanoplot : skip nanoplot
     --skip_porechop : skip porechop
"""
}


workflow AUTOCYCLER_WF {
     main:
     if (!params.outdir) { exit 1, "Missing parameter 'outdir'." }
     if (!params.input)  { exit 1, "Missing parameter 'input'."  }

     ch_versions = channel.empty()

     // 1. Validate sample sheet and dispatch samples
     ch_dispatched = DISPATCH().samplesheet

     // 2. Extract long reads for processing, map them to standard meta structure
     ch_raw_long_reads = ch_dispatched
          .map { row -> [ [ sample_id: row.sample_id ], row.long_reads ] }

     // 3. Run QC on long reads
     qc_out = LONGREAD_QC(ch_raw_long_reads)
     ch_versions = ch_versions.mix(qc_out.versions)

     ch_qc_reads = qc_out.long_reads_qc

     // 4. Calculate genome size for each QC'd sample
     ch_genome_size = GENOME_SIZE(ch_qc_reads)
     ch_versions = ch_versions.mix(GENOME_SIZE.out.versions)

     // 5. Subsample the reads
     ch_subsample_input = ch_qc_reads.join(ch_genome_size.size)
     ch_subsampled_reads = SUBSAMPLE(ch_subsample_input)
     ch_versions = ch_versions.mix(SUBSAMPLE.out.versions)

     // Flatten subsamples so each read subset can be processed independently
     ch_flat_subsamples = ch_subsampled_reads.reads.transpose()

     // 6. Assemble each subsampled file with multiple assemblers
     def assemblers = channel.fromList([
          'flye', 'plassembler', 'canu', 'metamdbg',
          'miniasm', 'myloasm', 'raven', 'necat', 'nextdenovo'
     ])

     ch_assembly_inputs = ch_flat_subsamples
          .map { meta, read_file -> [ meta.sample_id, read_file ] }
          .combine(assemblers)
          .combine(ch_genome_size.size.map { meta, size -> [meta.sample_id, size] }, by: 0)
          .map { sample_id, read_file, assembler, size ->
               def subsample_id = read_file.baseName.tokenize('_').last()
               def meta = [ id: "${sample_id}_subsample_${subsample_id}_${assembler}", sample_id: sample_id, assembler: assembler ]
               [ meta, read_file, size, params.read_type ]
          }

     ch_raw_assemblies = ASSEMBLE(ch_assembly_inputs)
     ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

     // 7. Apply clustering and consensus weights
     ch_weighted_assemblies = APPLY_WEIGHTS(ch_raw_assemblies.fasta)

     // Regroup weighted assemblies back by original sample ID
     ch_all_assemblies = ch_weighted_assemblies.fasta
          .map { meta, fasta -> [ meta.sample_id, fasta ] }
          .groupTuple()
          .map { sample_id, files -> [ [id: sample_id, sample_id: sample_id], files ] }

     // 8. Compress the input assemblies into a unitig graph
     ch_compressed = COMPRESS(ch_all_assemblies)
     ch_versions = ch_versions.mix(COMPRESS.out.versions)

     // 9. Cluster the input contigs into putative genomic sequences
     ch_clustered = CLUSTER(ch_compressed.graph)
     ch_versions = ch_versions.mix(CLUSTER.out.versions)

     // Extract individual QC-pass clusters (e.g., cluster_1, cluster_2) for parallel processing
     ch_qc_clusters = ch_clustered.qc_clusters
          .transpose()
          .map { meta, cluster_dir -> [ [id: "${meta.sample_id}_${cluster_dir.baseName}", sample_id: meta.sample_id], cluster_dir ] }

     // 10. Trim and resolve each QC-pass cluster
     ch_resolved_clusters = TRIM_AND_RESOLVE(ch_qc_clusters)
     ch_versions = ch_versions.mix(TRIM_AND_RESOLVE.out.versions)

     // Collect all resolved 5_final.gfa graphs to combine per sample
     ch_all_resolved = ch_resolved_clusters.gfa
          .map { meta, gfa -> [ meta.sample_id, gfa ] }
          .groupTuple()
          .map { sample_id, files -> [ [id: sample_id, sample_id: sample_id], files ] }

     // 11. Combine resolved clusters into a final assembly
     ch_combine_input = ch_all_resolved.join(ch_clustered.autocycler_out)
     ch_final_assembly = COMBINE(ch_combine_input)
     ch_versions = ch_versions.mix(COMBINE.out.versions)

     // Aggregate versions
     ch_versions = ch_versions.unique { version -> version.text }.collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")

     emit:
     nanoplot_raw_html   = qc_out.nanoplot_raw_html
     nanoplot_qc_html    = qc_out.nanoplot_qc_html
     porechop_log        = qc_out.porechop_log
     qc_long_reads       = ch_qc_reads
     final_assembly      = ch_final_assembly.fasta
     versions            = ch_versions
}

workflow {
    main:
    if (params.help) {
        helpMessage()
        exit(0, "")
    } else {
        AUTOCYCLER_WF()
    }

    publish:
    nanoplot_raw_html   = AUTOCYCLER_WF.out.nanoplot_raw_html
    nanoplot_qc_html    = AUTOCYCLER_WF.out.nanoplot_qc_html
    porechop_log        = AUTOCYCLER_WF.out.porechop_log
    qc_long_reads       = AUTOCYCLER_WF.out.qc_long_reads
    final_assembly      = AUTOCYCLER_WF.out.final_assembly
    versions            = AUTOCYCLER_WF.out.versions

}

output {
    qc_long_reads {
        path { meta, _file -> "sample_${meta.sample_id}/reads/post_qc" }
        mode "copy"
    }
    nanoplot_raw_html {
        path { meta, _file -> "sample_${meta.sample_id}/quality_assessment/nanoplot/raw" }
        mode "copy"
    }
    nanoplot_qc_html {
        path { meta, _file -> "sample_${meta.sample_id}/quality_assessment/nanoplot/post_qc" }
        mode "copy"
    }
    porechop_log {
        path { meta, _file -> "sample_${meta.sample_id}/quality_control/porechop" }
        mode "copy"
    }
    final_assembly {
        path { meta, _file -> "sample_${meta.sample_id}/autocycler_final" }
        mode "copy"
    }
    versions {
        path "pipeline_info/software_versions.yml"
        mode "copy"
    }
}
