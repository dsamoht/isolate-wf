process DNAAPLER {

    label 'process_low'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dnaapler:1.4.0--pyhdfd78af_0' :
        'biocontainers/dnaapler:1.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("dnaapler_out/*_reoriented.fasta"), emit: fasta
    path "versions.yml"                                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Only reorient circular contigs: list the linear ones for dnaapler to ignore
    grep '^>' ${assembly} | grep -v 'circular=true' | sed 's/^>//; s/ .*//' > linear_contigs.txt || true
    ignore_opt=""
    if [ -s linear_contigs.txt ]; then
        ignore_opt="--ignore linear_contigs.txt"
    fi

    dnaapler all \\
        ${args} \\
        -i ${assembly} \\
        -o dnaapler_out \\
        -p ${prefix} \\
        -t ${task.cpus} \\
        \$ignore_opt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaapler: \$(dnaapler --version | sed 's/.*version //')
    END_VERSIONS
    """
}
