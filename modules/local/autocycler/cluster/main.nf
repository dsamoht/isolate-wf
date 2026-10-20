process AUTOCYCLER_CLUSTER {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    tuple val(meta), path(autocycler_out)

    output:
    tuple val(meta), path("autocycler_out_work/clustering/qc_pass/cluster_*"), emit: qc_clusters, optional: true
    tuple val(meta), path("autocycler_out_work")                             , emit: autocycler_out
    path "versions.yml"                                                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cp -rL ${autocycler_out} autocycler_out_work

    autocycler cluster ${args} -a autocycler_out_work

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
