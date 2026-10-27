process SPLIT_TARGETS {
    tag "${meta.id}"
    label 'process_single'

    container "python:3.15.1-slim-trixie"

    input:
    tuple val(meta), path(bed)

    output:
    path("*_split.bed"), emit: bed
    tuple val("${task.process}"), val('split_targets.py'), val("v1.0.0"), topic: versions, emit: versions_extend_targets_py //TODO
    tuple val("${task.process}"), val('python'), eval("python -V | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    split_targets.py \\
        --input_targets_bed ${bed} \\
        --output_targets_bed ${prefix}_split.bed \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_split.bed
    """
}
