process CAT_FASTQ {

    label "small"

    tag "group_${meta.sample_id}"

    container params.coreutils_container

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    script:
    def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
    
    if (readList.size >= 1) {
    """
    cat ${readList.join(' ')} > ${meta.sample_id}.merged.fastq.gz
    """
    }
}
