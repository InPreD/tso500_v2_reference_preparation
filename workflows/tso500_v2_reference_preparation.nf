/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EXPAND_TARGETS                      } from '../subworkflows/local/expand_targets'
include { NORMAL_GATK_PREPROCESSING           } from '../subworkflows/local/normal_gatk_preprocessing'
include { paramsSummaryMap                    } from 'plugin/nf-schema'
include { softwareVersionsToYAML              } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TSO500_V2_REFERENCE_PREPARATION {

    take:
    ch_samplesheet        // channel: samplesheet read in from --input
    gatk4_resource_labels // string: semicolon-separated list of resource labels for gatk4 variantrecalibrator
    gatk4_resources       // string: path to directory containing gatk4 resource files for gatk4 variantrecalibrator
    minimum_depth         // integer: minimum depth for filter_hets_tsv
    panel_bed             // string: path to panel bed file
    reference             // string: path to reference file
    reference_dict        // string: path to reference file dict
    reference_fai         // string: path to reference file fai
    outdir                // string: path to output directory

    main:

    // general channels for versions and reference files
    def ch_versions = channel.empty()
    def ch_reference_simple = channel.fromPath(reference)
        .first()
    def ch_reference = ch_reference_simple
        .map{ fasta ->
            [ [ id: fasta.simpleName ], fasta ]
        }
    def ch_reference_dict_simple = channel.fromPath(reference_dict)
        .first()
    def ch_reference_dict = ch_reference_dict_simple
        .map{ dict ->
            [ [ id: dict.simpleName ], dict ]
        }
    def ch_reference_fai_simple = channel.fromPath(reference_fai)
        .first()
    def ch_reference_fai = ch_reference_fai_simple
        .map{ fai ->
            [ [ id: fai.simpleName ], fai ]
        }

    //
    // SUBWORKFLOW: Run Expand Targets
    //
    EXPAND_TARGETS(
        ch_reference_fai,
        ch_samplesheet,
        panel_bed
    )

    //
    // SUBWORKFLOW: Run Normal GATK Preprocessing
    //
    NORMAL_GATK_PREPROCESSING(
        EXPAND_TARGETS.out.bed,
        ch_reference,
        ch_reference_simple,
        ch_reference_dict,
        ch_reference_dict_simple,
        ch_reference_fai,
        ch_reference_fai_simple,
        ch_samplesheet,
        gatk4_resource_labels,
        gatk4_resources,
        minimum_depth
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        // remove consecutively repeated items
        .distinct()
        // divide items into two channels whether they are a Path or not
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        // map process to tool and version
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        // combine tool and version string according to process
        .groupTuple(by:0)
        // remove duplicates and sort tool and version strings, then format process: tool: version as YAML
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        // concatenate topic version strings
        .mix(topic_versions_string)
        // write versions to a single file
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
