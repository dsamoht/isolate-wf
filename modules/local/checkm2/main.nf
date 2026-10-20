process CHECKM2 {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/checkm2:1.1.0"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("checkm2_out/quality_report.tsv"), emit: report
    path "versions.yml"                                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Name the genome after the sample, which CheckM2 uses as its row name in the report
    mkdir genome
    ln -s ../${assembly} genome/${prefix}.fasta

    checkm2 predict \\
        ${args} \\
        --threads ${task.cpus} \\
        --input genome \\
        -x fasta \\
        --output-directory checkm2_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
