process CLUSTER {

    tag "${meta.sample_id}"

    label 'medium'

    container params.autocycler_container

    input:
    tuple val(meta), path(autocycler_out)

    output:
    tuple val(meta), path("autocycler_out_work/clustering/qc_pass/cluster_*"), emit: qc_clusters, optional: true
    tuple val(meta), path("autocycler_out_work")                             , emit: autocycler_out
    path "versions.yml"                                                      , emit: versions

    script:
    """
    cp -rL $autocycler_out autocycler_out_work

    autocycler cluster -a autocycler_out_work

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
