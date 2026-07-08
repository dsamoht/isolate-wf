process SUBSAMPLE {

    tag "${meta.sample_id}"
    
    label 'small'

    container params.autocycler_container

    input:
    tuple val(meta), path(reads), val(genome_size)

    output:
    tuple val(meta), path("subsampled_reads/*.fastq"), emit: reads
    path "versions.yml"                              , emit: versions

    script:
    """
    autocycler subsample \\
        --reads $reads \\
        --out_dir subsampled_reads \\
        --genome_size $genome_size

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
