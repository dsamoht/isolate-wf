process SKANI_TRIANGLE {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skani:0.3.2--h79ce301_0' :
        'biocontainers/skani:0.3.2--h79ce301_0' }"

    input:
    // Every final assembly of the run, with the sample ids in the same order. Each assembly
    // is staged in a directory of its own: without Dnaapler they are all consensus_assembly.fasta
    tuple val(meta), val(ids), path(genomes, stageAs: 'input*/*')

    output:
    tuple val(meta), path("skani_triangle.tsv")  , emit: ani
    tuple val(meta), path("duplicate_isolates.tsv"), emit: duplicates
    path "versions.yml"                          , emit: versions

    script:
    def args    = task.ext.args ?: ''
    // Two isolates are flagged as the same strain above both thresholds; the aligned
    // fraction must pass in both directions
    def min_ani = task.ext.min_ani != null ? task.ext.min_ani : 99.9
    def min_af  = task.ext.min_af  != null ? task.ext.min_af  : 95
    def links   = [ids, genomes].transpose().collect { id, genome -> "ln -s ../${genome} genomes/${id}.fasta" }.join('\n    ')
    """
    mkdir -p genomes
    ${links}

    # Every pair, as an edge list: Ref_file, Query_file, ANI, Align_fraction_ref,
    # Align_fraction_query, Ref_name, Query_name
    skani triangle \\
        ${args} \\
        -E \\
        -t ${task.cpus} \\
        genomes/*.fasta \\
        > skani_triangle.tsv

    printf 'sample_1\\tsample_2\\tani\\taf_1\\taf_2\\n' > duplicate_isolates.tsv
    awk -F '\\t' -v ani="${min_ani}" -v af="${min_af}" '
        NR > 1 && \$3 >= ani && \$4 >= af && \$5 >= af {
            n1 = \$1; sub(/^.*\\//, "", n1); sub(/\\.fasta\$/, "", n1)
            n2 = \$2; sub(/^.*\\//, "", n2); sub(/\\.fasta\$/, "", n2)
            print n1 "\\t" n2 "\\t" \$3 "\\t" \$4 "\\t" \$5
        }' skani_triangle.tsv >> duplicate_isolates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/skani //')
    END_VERSIONS
    """
}
