process TRIM_AND_RESOLVE {

    tag "${meta.sample_id}"
    
    label 'process_medium'

    container params.autocycler_container

    input:
    tuple val(meta), path(cluster_dir)

    output:
    tuple val(meta), path("${cluster_dir.name}/5_final.gfa"), emit: gfa
    path "versions.yml"                                     , emit: versions

    script:
    """
    # 1. Rename the Nextflow symlink temporarily to break the copy loop
    mv ${cluster_dir} ${cluster_dir}_symlink
    
    # 2. Copy the actual contents back into a real, mutable directory with the original name
    cp -rL ${cluster_dir}_symlink ${cluster_dir.name}
    
    # 3. Run Autocycler commands
    autocycler trim -c ${cluster_dir.name}
    autocycler resolve -c ${cluster_dir.name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
