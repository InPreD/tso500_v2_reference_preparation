process FAI_TO_CHR_SIZES {
    tag "${meta.id}"
    label 'process_single'

    container "ubuntu:24.04"

    input:
    tuple val(meta), path(fai)

    output:
    path("*_chr_sizes.tsv"), emit: tsv
    path("*_chr_sizes_sorted.tsv"), emit: sorted_tsv
    tuple val("${task.process}"), val('awk'), eval("awk --version | head -n1 | cut -d',' -f1 | cut -d' ' -f3"), topic: versions, emit: versions_sort

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk \\
        '{print \$1"\\t"\$2}' \\
        ${fai} \\
        > ${prefix}_chr_sizes.tsv
    sort \\
        -k1,1 \\
        -k2,2n \\
        ${prefix}_chr_sizes.tsv \\
        > ${prefix}_chr_sizes_sorted.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chr_sizes.tsv
    touch ${prefix}_chr_sizes_sorted.tsv
    """
}
