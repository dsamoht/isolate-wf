process AUTOCYCLER_COMPRESS {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("autocycler_out"), emit: graph
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p assemblies_dir

    # Copy files resolving symlinks to prevent naming collisions and read-only issues
    for file in ${assemblies}; do
        cp -L "\$file" assemblies_dir/
    done

    autocycler compress \\
        ${args} \\
        --threads ${task.cpus} \\
        -i assemblies_dir \\
        -a autocycler_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
