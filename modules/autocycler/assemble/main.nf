process ASSEMBLE {

    tag "${meta.sample_id}"

    label 'medium'

    container params.autocycler_container

    input:
    tuple val(meta), path(reads), val(genome_size), val(read_type)

    output:
    tuple val(meta), path("${meta.assembler}_*.fasta"), emit: fasta, optional: true
    path "versions.yml"                               , emit: versions

    script:
    def sub_id = reads.baseName.tokenize('_').last()
    """
    autocycler helper ${meta.assembler} \\
        --reads $reads \\
        --out_prefix ${meta.assembler}_${sub_id} \\
        --threads $task.cpus \\
        --genome_size $genome_size \\
        --read_type $read_type \\
        --min_depth_rel 0.1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
