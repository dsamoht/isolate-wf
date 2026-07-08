process COMBINE {

    tag "${meta.sample_id}"

    label 'process_medium'

    container params.autocycler_container

    input:
    tuple val(meta), path(resolved_gfas), path(autocycler_out)

    output:
    tuple val(meta), path("autocycler_out_final/consensus_assembly.fasta"), emit: fasta, optional: true
    path "versions.yml"                                                   , emit: versions

    script:
    """
    # Copy the core directory structure
    cp -rL $autocycler_out autocycler_out_final
    
    autocycler combine \\
        -a autocycler_out_final \\
        -i $resolved_gfas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
