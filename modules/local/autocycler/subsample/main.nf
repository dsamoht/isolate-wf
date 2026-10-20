process AUTOCYCLER_SUBSAMPLE {

    label 'process_low'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    tuple val(meta), path(reads), val(genome_size)

    output:
    tuple val(meta), path("subsampled_reads/*.fastq"), emit: reads
    path "versions.yml"                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    autocycler subsample \\
        ${args} \\
        --reads ${reads} \\
        --out_dir subsampled_reads \\
        --genome_size ${genome_size}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
