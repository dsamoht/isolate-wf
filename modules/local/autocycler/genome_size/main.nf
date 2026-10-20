process AUTOCYCLER_GENOME_SIZE {

    label 'process_low'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), env('GENOME_SIZE'), emit: size
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    GENOME_SIZE=\$(autocycler helper genome_size ${args} --reads ${reads} --threads ${task.cpus})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
