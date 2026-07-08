process COMPRESS {

    tag "${meta.sample_id}"
    
    label 'process_medium'

    container params.autocycler_container

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("autocycler_out"), emit: graph
    path "versions.yml"                    , emit: versions

    script:
    """
    mkdir -p assemblies_dir
    
    # Copy files resolving symlinks to prevent naming collisions and read-only issues
    for file in $assemblies; do
        cp -L "\$file" assemblies_dir/
    done
    
    autocycler compress \\
        -i assemblies_dir \\
        -a autocycler_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
