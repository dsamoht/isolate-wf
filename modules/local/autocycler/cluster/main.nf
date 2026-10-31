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
    // Emitted for every sample, also when no cluster passes QC, so that case can be inspected
    tuple val(meta), path("clustering.{tsv,newick,yaml}")                    , emit: clustering
    // The whole Autocycler directory, for the samples AUTOCYCLER_COMBINE never sees: with
    // its input_assemblies.gfa, `autocycler cluster --manual` can be re-run on it
    tuple val(meta), path("*_autocycler_out.tar.gz")                         , emit: archive
    path "versions.yml"                                                      , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -rL ${autocycler_out} autocycler_out_work

    autocycler cluster ${args} -a autocycler_out_work

    # Copied out of the Autocycler directory so they publish without the nesting; they stay
    # in it too, for AUTOCYCLER_COMBINE and the archive
    cp autocycler_out_work/clustering/clustering.{tsv,newick,yaml} .

    # Archived as autocycler_out/, the name Autocycler gives it, through a symlink that
    # `tar -h` follows: GNU and BSD tar both take it, unlike GNU's --transform
    mkdir archive
    ln -s ../autocycler_out_work archive/autocycler_out
    tar -chzf ${prefix}_autocycler_out.tar.gz -C archive autocycler_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
