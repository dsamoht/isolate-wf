process AUTOCYCLER_COMPRESS {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/autocycler:0.8.0"

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("autocycler_out")     , emit: graph
    // The contig count of every input assembly and whether it was kept, excluded or empty
    tuple val(meta), path("assembly_filter.tsv"), emit: filter
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    // The filter's floor follows `--max_contigs` when it is given in ext.args, so that the
    // filter never drops an assembly that Autocycler itself would have accepted
    def max_contigs_match = args =~ /--max_contigs[=\s]+(\d+)/
    def max_contigs       = max_contigs_match.find() ? max_contigs_match.group(1) as int : 25
    // Assemblies with more than max(min_limit, outlier_factor x median) contigs are excluded
    def outlier_factor = task.ext.outlier_factor    != null ? task.ext.outlier_factor    : 5
    def min_limit      = task.ext.outlier_min_limit != null ? task.ext.outlier_min_limit : max_contigs
    // Fail if fewer than this fraction of the non-empty assemblies are left after exclusion
    def min_kept       = task.ext.min_kept_fraction != null ? task.ext.min_kept_fraction : 0.5
    """
    mkdir -p assemblies_dir excluded_assemblies

    # Copy files resolving symlinks to prevent naming collisions and read-only issues
    for file in ${assemblies}; do
        cp -L "\$file" assemblies_dir/
    done

    # Exclude outlier assemblies, e.g. a metagenome assembler that keeps hundreds of
    # low-coverage contaminant contigs, so they do not push the mean contig count over
    # Autocycler's --max_contigs threshold for an otherwise clean isolate. Empty assemblies
    # (an assembler that ran but found nothing, as Plassembler does for a sample without
    # plasmids) are set aside too, as `autocycler compress` stops on an empty file, and left
    # out of the median and the counts, where they would only make the filter stricter.
    for f in assemblies_dir/*.fasta; do
        echo "\$(grep -c '>' "\$f" || true) \$f"
    done | sort -n > contig_counts.txt
    awk '\$1 > 0' contig_counts.txt > nonempty_counts.txt
    n_total=\$(awk 'END { print NR }' nonempty_counts.txt)
    if [ "\$n_total" -eq 0 ]; then
        echo "Error: none of the input assemblies has any contig" >&2
        exit 1
    fi
    median=\$(awk -v n="\$n_total" 'NR == int((n + 1) / 2) { print \$1 }' nonempty_counts.txt)
    limit=\$(awk -v m="\$median" -v f="${outlier_factor}" -v l="${min_limit}" 'BEGIN { x = int(m * f); print (x > l ? x : l) }')
    echo "Contig count per assembly: median \$median over \$n_total non-empty assemblies, exclusion limit \$limit" >&2

    printf 'assembly\\tcontigs\\tlimit\\tstatus\\n' > assembly_filter.tsv
    n_kept=0
    while read -r count f; do
        if [ "\$count" -eq 0 ]; then
            status=empty
            mv "\$f" excluded_assemblies/
        elif [ "\$count" -gt "\$limit" ]; then
            status=excluded
            echo "Excluding \$(basename "\$f"): \$count contigs" >&2
            mv "\$f" excluded_assemblies/
        else
            status=kept
            n_kept=\$((n_kept + 1))
        fi
        printf '%s\\t%s\\t%s\\t%s\\n' "\$(basename "\$f" .fasta)" "\$count" "\$limit" "\$status" >> assembly_filter.tsv
    done < contig_counts.txt

    if ! awk -v k="\$n_kept" -v t="\$n_total" -v m="${min_kept}" 'BEGIN { exit !(k >= m * t) }'; then
        echo "Error: only \$n_kept of \$n_total non-empty assemblies left after excluding outliers (minimum fraction ${min_kept})" >&2
        exit 1
    fi

    autocycler compress \\
        ${args} \\
        --threads ${task.cpus} \\
        -i assemblies_dir \\
        -a autocycler_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | sed 's/Autocycler //g')
    END_VERSIONS
    """
}
