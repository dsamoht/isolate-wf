process GENOMAD {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.13.0--pyhdfd78af_0' :
        'biocontainers/genomad:1.13.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)
    path genomad_db

    output:
    tuple val(meta), path("genomad_out"), emit: results
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # geNomad names its output files after the input file, so give it the sample name
    ln -s ${assembly} ${prefix}.fasta

    genomad end-to-end \\
        ${args} \\
        --cleanup \\
        --threads ${task.cpus} \\
        ${prefix}.fasta \\
        genomad_out \\
        ${genomad_db}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(genomad --version | sed 's/.*version //')
    END_VERSIONS
    """
}
