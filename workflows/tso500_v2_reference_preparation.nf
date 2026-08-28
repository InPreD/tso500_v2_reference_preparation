/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BEDTOOLS_BAMTOBED      } from '../modules/nf-core/bedtools/bamtobed'
include { BEDTOOLS_COVERAGE      } from '../modules/nf-core/bedtools/coverage'
include { BEDTOOLS_MERGE         } from '../modules/nf-core/bedtools/merge'
include { BEDTOOLS_SLOP          } from '../modules/nf-core/bedtools/slop'
include { BEDTOOLS_SUBTRACT      } from '../modules/nf-core/bedtools/subtract'
include { FAI_TO_CHR_SIZES       } from '../modules/local/fai_to_chr_sizes'
include { SORT as SORT_BED       } from '../modules/local/sort'
include { SORT as SORT_SUBTRACT  } from '../modules/local/sort'
include { SORT as SORT_BAMTOBED  } from '../modules/local/sort'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TSO500_V2_REFERENCE_PREPARATION {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    panel_bed
    reference_fai
    outdir

    main:

    def ch_versions = channel.empty()
    //
    // MODULE: Run fai to chromosome sizes
    //
    ch_fai_to_chr_sizes_input = channel.fromPath(reference_fai).map{ fai -> [ [ id: fai.simpleName ], fai ] }
    FAI_TO_CHR_SIZES(ch_fai_to_chr_sizes_input)

    //
    // MODULE: Run sort
    //
    ch_sort_input = channel.fromPath(panel_bed).map{ bed -> [ [ id: bed.simpleName ], bed ] }
    SORT_BED(ch_sort_input)

    //
    // MODULE: Run bedtools/slop
    //
    ch_bedtools_slop_input = SORT_BED.out.bed.map{ meta, bed -> [ [ id: meta.id + '_slop' ], bed ] }
    BEDTOOLS_SLOP(ch_bedtools_slop_input, FAI_TO_CHR_SIZES.out.tsv)

    //
    // MODULE: Run bedtools/subtract
    //
    ch_bedtools_subtract_input = BEDTOOLS_SLOP.out.bed.map{ meta, bed -> [ [ id: meta.id + '_complement' ], bed, file(panel_bed) ] }
    BEDTOOLS_SUBTRACT(ch_bedtools_subtract_input)

    //
    // MODULE: Run sort
    //
    SORT_SUBTRACT(BEDTOOLS_SUBTRACT.out.bed)

    //
    // MODULE: Run bedtools/merge
    //
    ch_bedtools_merge_input = SORT_SUBTRACT.out.bed.map{ meta, bed -> [ [ id: meta.id + '_merged' ], bed ] }
    BEDTOOLS_MERGE(ch_bedtools_merge_input)

    //
    // MODULE: Run bedtools/bamtobed
    //
    BEDTOOLS_BAMTOBED(ch_samplesheet)

    //
    // MODULE: Run sort
    //
    SORT_BAMTOBED(BEDTOOLS_BAMTOBED.out.bed)

    //
    // MODULE: Run bedtools/coverage
    //
    ch_bedtools_coverage_input = SORT_BAMTOBED.out.bed.combine(BEDTOOLS_MERGE.out.bed).map{ meta1, bed1, meta2, bed2 -> [ [ id: meta1.id + '_' + meta2.id ], bed2, bed1 ] }
    BEDTOOLS_COVERAGE(ch_bedtools_coverage_input, FAI_TO_CHR_SIZES.out.sorted_tsv)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'tso500_v2_reference_preparation_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
