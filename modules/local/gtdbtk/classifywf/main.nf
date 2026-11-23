process GTDBTK_CLASSIFYWF {

    label 'process_high'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.6.1--pyh1f0d9b5_2' :
        'biocontainers/gtdbtk:2.6.1--pyh1f0d9b5_2' }"

    input:
    // Every final assembly of the run, with the sample ids in the same order. Each assembly
    // is staged in a directory of its own: without Dnaapler they are all consensus_assembly.fasta
    tuple val(meta), val(ids), path(genomes, stageAs: 'input*/*')
    path gtdbtk_db

    output:
    tuple val(meta), path("gtdbtk_out")           , emit: results
    tuple val(meta), path("gtdbtk.*.summary.tsv") , emit: summary
    path "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    // GTDB-Tk names each genome after its file: one file per sample id
    def links = [ids, genomes].transpose().collect { id, genome -> "ln -s ../${genome} genomes/${id}.fasta" }.join('\n    ')
    """
    mkdir -p genomes
    ${links}

    export GTDBTK_DATA_PATH=\$(readlink -f ${gtdbtk_db})

    gtdbtk classify_wf \\
        ${args} \\
        --genome_dir genomes \\
        --extension fasta \\
        --out_dir gtdbtk_out \\
        --cpus ${task.cpus} \\
        --pplacer_cpus ${Math.min(task.cpus as int, 8)}

    # One table per domain found among the genomes: bac120, ar53 or both
    cp gtdbtk_out/gtdbtk.*.summary.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version 2>&1 | sed -n 's/.*version \\([0-9.]*\\).*/\\1/p' | head -n 1)
    END_VERSIONS
    """
}
