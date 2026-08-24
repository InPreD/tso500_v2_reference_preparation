process SORT {
    tag "${meta.id}"
    label 'process_single'

    container "ubuntu:24.04"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('sort'), eval("sort --version | grep '^sort' | sed -e 's/sort (GNU coreutils) //g'"), topic: versions, emit: versions_sort

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sort \\
        ${args} \\
        ${bed} \\
        > ${prefix}_sort.bed
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    touch ${prefix}_sort.bed
    """
}
