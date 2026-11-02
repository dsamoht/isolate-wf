process SAMPLE_SUMMARY {

    label 'process_single'

    tag "$meta.id"

    // Only the Python standard library is needed: the Autocycler image, which every run
    // pulls anyway, has Python 3
    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    // One directory per step under inputs/; a step that was skipped or never reached for
    // this sample comes as an empty list and leaves its directory out
    tuple val(meta), val(info),
        path(nanoplot_raw, stageAs: 'inputs/nanoplot_raw/*'),
        path(porechop    , stageAs: 'inputs/porechop/*'),
        path(chopper     , stageAs: 'inputs/chopper/*'),
        path(nanoplot_qc , stageAs: 'inputs/nanoplot_qc/*'),
        path(subsample   , stageAs: 'inputs/subsample/*'),
        path(assemblies  , stageAs: 'inputs/assemblies/*'),
        path(clustering  , stageAs: 'inputs/clustering/*'),
        path(metrics     , stageAs: 'inputs/metrics/*'),
        path(consensus   , stageAs: 'inputs/consensus/*'),
        path(dnaapler    , stageAs: 'inputs/dnaapler/*'),
        path(checkm2     , stageAs: 'inputs/checkm2/*'),
        path(genomad     , stageAs: 'inputs/genomad/*'),
        path(versions)

    output:
    tuple val(meta), path("isolate-wf_sample_*_summary.log"), emit: log

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // One `key<TAB>value` line per entry, each single-quoted for the shell. Written with
    // printf rather than a heredoc: the interpolated lines would carry no indentation, and
    // Nextflow would then no longer strip it from the heredoc's closing line
    def info_lines = info.collect { key, value ->
        def text = value instanceof List ? value.join(',') : (value == null ? '' : value.toString())
        "'" + "${key}\t${text}".replace("'", "'\\''") + "'"
    }.join(' ')
    """
    printf '%s\\n' ${info_lines} > info.tsv

    isolate_wf_summary.py \\
        --info info.tsv \\
        --versions ${versions} \\
        --inputs inputs \\
        --output isolate-wf_sample_${prefix}_summary.log
    """
}
