process FILTER_HETS_TSV {
    tag "${meta.id}"
    label 'process_single'

    container "ubuntu:24.04"

    input:
    tuple val(meta), path(tsv)
    val(minimum_depth)

    output:
    path("*_filtered.hets.tsv"), emit: tsv
    tuple val("${task.process}"), val('awk'), eval("awk --version | head -n1 | cut -d',' -f1 | cut -d' ' -f3"), topic: versions, emit: versions_awk
    tuple val("${task.process}"), val('grep'), eval("grep --version | head -n1 | cut -d' ' -f4"), topic: versions, emit: versions_grep

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep "^@\\|CONTIG" ${tsv} \\
        > ${prefix}_filtered.hets.tsv
    grep -v "^@\\|CONTIG" ${tsv} \\
        | awk \\
            -v minimum=${minimum_depth} \\
            '{if ((\$3+\$4) >= minimum) print \$0}'
            >> ${prefix}_filtered.hets.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered.hets.tsv
    """
}
