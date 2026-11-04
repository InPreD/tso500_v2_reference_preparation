process POSITION_COUNT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(bed)
    path tsv

    output:
    tuple val(meta), path("*_position_count.tsv"), emit: tsv
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "file\\tposition_count_in_millions" >> ${prefix}_position_count.tsv
    for FILE in ${bed}; do
        COUNT=`bedtools \\
            genomecov \\
            ${args} \\
            -i \$FILE \\
            -g ${tsv} \\
            | wc -l \\
            | awk '{print \$1/1000000}'`;
        echo -e "\$FILE\\t\$COUNT" >> ${prefix}_position_count.tsv;
    done
    
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_position_count.tsv
    """
}
