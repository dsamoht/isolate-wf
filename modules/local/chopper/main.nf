process CHOPPER {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.9.0--hdcf5f25_0' :
        'biocontainers/chopper:0.9.0--hdcf5f25_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.choppered.fastq.gz"), emit: fastq
    // chopper reports how many reads it kept on stderr
    tuple val(meta), path("*.chopper.log")       , emit: log
    path "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # chopper reads gzipped and plain FastQ alike when given the file itself
    # Its stderr goes to the log, then on to the task's own stderr, failed or not, so an
    # error still shows in .command.err
    chopper ${args} --threads ${task.cpus} -i ${reads} 2> ${prefix}.chopper.log \
        | gzip > ${prefix}.choppered.fastq.gz \
        || { cat ${prefix}.chopper.log >&2; exit 1; }
    cat ${prefix}.chopper.log >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
