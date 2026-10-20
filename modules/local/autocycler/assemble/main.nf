process AUTOCYCLER_ASSEMBLE {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    tuple val(meta), path(reads), val(genome_size), val(read_type)

    output:
    // An assembler that fails leaves no FASTA: emit an empty list rather than nothing, so
    // the per-sample regroup downstream still receives one item per assembly task
    tuple val(meta), path("${meta.assembler}_*.fasta", arity: '0..*'), emit: fasta
    path "versions.yml"                                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    // The assembler and subsample name the file, which `autocycler compress` reads back
    def sub_id = reads.baseName.tokenize('_').last()
    def prefix = "${meta.assembler}_${sub_id}"
    """
    autocycler helper ${meta.assembler} \\
        ${args} \\
        --reads ${reads} \\
        --out_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --genome_size ${genome_size} \\
        --read_type ${read_type}

    # Clustering and consensus weights, as recommended in the Autocycler docs
    if [[ -s ${prefix}.fasta ]]; then
        if [[ "${meta.assembler}" == "plassembler" ]]; then
            sed -i 's/circular=True/circular=True Autocycler_cluster_weight=3/' ${prefix}.fasta
        elif [[ "${meta.assembler}" == "canu" || "${meta.assembler}" == "flye" ]]; then
            sed -i 's/^>.*\$/& Autocycler_consensus_weight=2/' ${prefix}.fasta
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
