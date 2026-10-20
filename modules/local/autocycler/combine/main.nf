process AUTOCYCLER_COMBINE {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    tuple val(meta), path(resolved_gfas), path(autocycler_out)

    output:
    tuple val(meta), path("consensus_assembly.fasta"), emit: fasta, optional: true
    path "versions.yml"                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Copy the core directory structure
    cp -rL ${autocycler_out} autocycler_out_final

    autocycler combine \\
        ${args} \\
        --threads ${task.cpus} \\
        -a autocycler_out_final \\
        -i ${resolved_gfas}

    # Lift the consensus out of the Autocycler directory so it publishes without the nesting
    if [ -f autocycler_out_final/consensus_assembly.fasta ]; then
        mv autocycler_out_final/consensus_assembly.fasta .
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
