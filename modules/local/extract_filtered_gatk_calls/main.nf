process EXTRACT_FILTERED_GATK_CALLS {
    tag "${meta.id}"
    label 'process_single'

    container "inpred/python:3.15.2-slim-trixie"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('extract_filtered_gatk_calls.py'), val("${workflow.manifest.version}"), topic: versions, emit: versions_extract_filtered_gatk_calls_py
    tuple val("${task.process}"), val('python'), eval("python -V | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_filtered_gatk_calls.py \\
        --input_vcf ${vcf} \\
        --output_tsv ${prefix}.bed \\
        --output_vcf ${prefix}.vcf \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed ${prefix}.vcf
    """
}
