process EXTEND_TARGETS {
    tag "${meta1.id}"
    label 'process_single'

    container "python:3.15.2-slim-trixie"

    input:
    tuple val(meta1), path(bed)
    tuple val(meta2), path(tsv)

    output:
    path("*_grown.bed"), emit: bed
    tuple val("${task.process}"), val('extend_targets.py'), val("${workflow.manifest.version}"), topic: versions, emit: versions_extend_targets_py
    tuple val("${task.process}"), val('python'), eval("python -V | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta1.id}"
    def coverage_tsv = tsv.join(' --coverage_tsv ')
    """
    extend_targets.py \\
        --input_targets_bed ${bed} \\
        --coverage_tsv ${coverage_tsv} \\
        --male_sample_count ${meta2.male} \\
        --female_sample_count ${meta2.female} \\
        --output_targets_bed ${prefix}_grown.bed \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta1.id}"
    """
    touch ${prefix}_grown.bed
    """
}
