process AUTOCYCLER_COMBINE {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    // Cluster directory names (cluster_001, ...) are unique within a sample, so the resolved
    // clusters can share one staging directory
    tuple val(meta), path(resolved_clusters, stageAs: "resolved_clusters/*"), path(autocycler_out), path(subsample_yaml)

    output:
    tuple val(meta), path("consensus_assembly.fasta")     , emit: fasta  , optional: true
    tuple val(meta), path("consensus_assembly.gfa")       , emit: gfa    , optional: true
    tuple val(meta), path("*_autocycler_metrics.tsv")     , emit: metrics
    path "versions.yml"                                   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Rebuild a complete Autocycler directory, laid out as a run of autocycler_full.sh
    # would leave it: `autocycler table` searches it for the YAML metrics of every step
    cp -rL ${autocycler_out} autocycler_out_final
    for cluster in resolved_clusters/*; do
        rm -rf autocycler_out_final/clustering/qc_pass/\$(basename \$cluster)
        cp -rL \$cluster autocycler_out_final/clustering/qc_pass/
    done
    cp -L ${subsample_yaml} autocycler_out_final/

    autocycler combine \\
        ${args} \\
        --threads ${task.cpus} \\
        -a autocycler_out_final \\
        -i autocycler_out_final/clustering/qc_pass/cluster_*/5_final.gfa

    # Header line, then this sample's line
    {
        autocycler table
        autocycler table -a autocycler_out_final -n ${prefix}
    } > ${prefix}_autocycler_metrics.tsv

    # Lift the consensus out of the Autocycler directory so it publishes without the nesting
    for ext in fasta gfa; do
        if [ -f autocycler_out_final/consensus_assembly.\$ext ]; then
            mv autocycler_out_final/consensus_assembly.\$ext .
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
