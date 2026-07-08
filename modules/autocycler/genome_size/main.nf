process GENOME_SIZE {

    tag "${meta.sample_id}"

    label 'small'
    
    container params.autocycler_container

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), env(GENOME_SIZE), emit: size
    path "versions.yml"              , emit: versions

    script:
    """
    GENOME_SIZE=\$(autocycler helper genome_size --reads $reads --threads $task.cpus)
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
