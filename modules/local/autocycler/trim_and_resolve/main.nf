process AUTOCYCLER_TRIM_AND_RESOLVE {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    tuple val(meta), path(cluster_dir)

    output:
    tuple val(meta), path("${cluster_dir.name}/5_final.gfa"), emit: gfa
    path "versions.yml"                                     , emit: versions

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    # 1. Rename the Nextflow symlink temporarily to break the copy loop
    mv ${cluster_dir} ${cluster_dir}_symlink

    # 2. Copy the actual contents back into a real, mutable directory with the original name
    cp -rL ${cluster_dir}_symlink ${cluster_dir.name}

    # 3. Run Autocycler commands
    autocycler trim ${args} --threads ${task.cpus} -c ${cluster_dir.name}
    autocycler resolve ${args2} -c ${cluster_dir.name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
