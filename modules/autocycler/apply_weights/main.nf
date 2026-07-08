process APPLY_WEIGHTS {

    tag "${meta.sample_id}"
    
    label 'small'

    container params.coreutils_container

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("weighted_${fasta.name}"), emit: fasta

    script:
    """
    cp -L $fasta weighted_${fasta.name}

    if [[ "${meta.assembler}" == "plassembler" ]]; then
        sed -i 's/circular=True/circular=True Autocycler_cluster_weight=3/' weighted_${fasta.name}
    elif [[ "${meta.assembler}" == "canu" || "${meta.assembler}" == "flye" ]]; then
        sed -i 's/^>.*\$/& Autocycler_consensus_weight=2/' weighted_${fasta.name}
    fi
    """
}
