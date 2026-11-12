/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_NORM                             } from '../modules/nf-core/bcftools/norm'
include { BEDTOOLS_BAMTOBED                         } from '../modules/nf-core/bedtools/bamtobed'
include { BEDTOOLS_COVERAGE                         } from '../modules/nf-core/bedtools/coverage'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_SUBTRACT } from '../modules/nf-core/bedtools/merge'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_EXTEND   } from '../modules/nf-core/bedtools/merge'
include { BEDTOOLS_SLOP                             } from '../modules/nf-core/bedtools/slop'
include { BEDTOOLS_SUBTRACT                         } from '../modules/nf-core/bedtools/subtract'
include { EXTEND_TARGETS                            } from '../modules/local/extend_targets'
include { EXTRACT_FILTERED_GATK_CALLS               } from '../modules/local/extract_filtered_gatk_calls'
include { FAI_TO_CHR_SIZES                          } from '../modules/local/fai_to_chr_sizes'
include { GATK4_ADDORREPLACEREADGROUPS              } from '../modules/nf-core/gatk4/addorreplacereadgroups'
include { GATK4_APPLYVQSR                           } from '../modules/nf-core/gatk4/applyvqsr'
include { GATK4_COLLECTREADCOUNTS                   } from '../modules/nf-core/gatk4/collectreadcounts'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS       } from '../modules/nf-core/gatk4/createreadcountpanelofnormals'
include { GATK4_DENOISEREADCOUNTS                   } from '../modules/nf-core/gatk4/denoisereadcounts'
include { GATK4_GENOMICSDBIMPORT                    } from '../modules/nf-core/gatk4/genomicsdbimport'
include { GATK4_GENOTYPEGVCFS                       } from '../modules/nf-core/gatk4/genotypegvcfs'
include { GATK4_HAPLOTYPECALLER                     } from '../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_PREPROCESSINTERVALS                 } from '../modules/nf-core/gatk4/preprocessintervals'
include { GATK4_VARIANTRECALIBRATOR                 } from '../modules/nf-core/gatk4/variantrecalibrator'
include { POSITION_COUNT                            } from '../modules/local/position_count'
include { SORT as SORT_BED                          } from '../modules/local/sort'
include { SORT as SORT_SUBTRACT                     } from '../modules/local/sort'
include { SORT as SORT_BAMTOBED                     } from '../modules/local/sort'
include { SPLIT_TARGETS                             } from '../modules/local/split_targets'
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TSO500_V2_REFERENCE_PREPARATION {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    panel_bed
    gatk4_resource_labels
    gatk4_resources
    reference
    reference_dict
    reference_fai
    outdir

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
    // MODULE: Run fai to chromosome sizes
    //
    FAI_TO_CHR_SIZES(ch_reference_fai)

    //
    // MODULE: Run sort
    //
    def ch_sort_input = channel.fromPath(panel_bed)
        // add meta object to bed file
        .map{ bed ->
            [ [ id: bed.simpleName ], bed ]
        }
    SORT_BED(ch_sort_input)

    //
    // MODULE: Run bedtools/slop
    //
    BEDTOOLS_SLOP(SORT_BED.out.bed, FAI_TO_CHR_SIZES.out.sorted_tsv.first())
    def ch_position_count_input = SORT_BED.out.bed
        // add bedtools slop output bed
        .combine(
            BEDTOOLS_SLOP.out.bed
                // remove meta object from bed file
                .map{ meta, bed ->
                    bed
                }
        )
        // combine bed files into single list
        .map{ meta, bed1, bed2 ->
            [ meta, [ bed1, bed2 ] ]
        }

    //
    // MODULE: Run bedtools/subtract
    //
    def ch_bedtools_subtract_input = BEDTOOLS_SLOP.out.bed
        // add panel bed to bedtools slop output bed
        .map{ meta, bed ->
            [ meta, bed, file(panel_bed) ]
        }
    BEDTOOLS_SUBTRACT(ch_bedtools_subtract_input)
    ch_position_count_input = ch_position_count_input
        // add bedtools subtract output bed
        .combine(
            BEDTOOLS_SUBTRACT.out.bed
                // remove meta object from bed file
                .map{ meta, bed ->
                    bed
                }
        )
        // combine bed files into single list
        .map{ meta, bed1, bed2 ->
            [ meta, bed1 + [ bed2 ] ]
        }

    //
    // MODULE: Run sort
    //
    SORT_SUBTRACT(BEDTOOLS_SUBTRACT.out.bed)

    //
    // MODULE: Run bedtools/merge
    //
    BEDTOOLS_MERGE_SUBTRACT(SORT_SUBTRACT.out.bed)
    ch_position_count_input = ch_position_count_input
        // add bedtools merge output bed
        .combine(
            BEDTOOLS_MERGE_SUBTRACT.out.bed
                // remove meta object from bed file
                .map{ meta, bed ->
                    bed
                }
        )
        // combine bed files into single list
        .map{ meta, bed1, bed2 ->
            [ meta, bed1 + [ bed2 ] ]
        }

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
    def ch_bedtools_coverage_input = SORT_BAMTOBED.out.bed
        // add bedtools merge output bed
        .combine(BEDTOOLS_MERGE_SUBTRACT.out.bed)
        // create combined meta object containing sex
        .map{ meta1, bed1, meta2, bed2 ->
            [ [ id: meta1.id + '_' + meta2.id, sex: meta1.sex.toLowerCase() ], bed2, bed1 ]
        }
    BEDTOOLS_COVERAGE(ch_bedtools_coverage_input, FAI_TO_CHR_SIZES.out.sorted_tsv.first())

    //
    // MODULE: Run extend_targets
    //
    def ch_extend_targets_input = BEDTOOLS_COVERAGE.out.bed
        // wrap bed file in list
        .map{ it ->
            [ it ]
        }
        // combine all bed files in channel into single list
        .collect()
        // count sex occurences and combine all bed files
        .map{ it ->
            def female = 0
            def male = 0
            def tsv = []
            it.each {
                meta, bed -> 
                if (meta.sex == 'f') {
                    female = female + 1
                } else {
                    male = male + 1
                }
                tsv.add(bed)
            }
            return [ [ id: 'normal', female: female, male: male ], tsv ]
        }
    EXTEND_TARGETS(ch_sort_input, ch_extend_targets_input)
    ch_position_count_input = ch_position_count_input
        // add extend targets output bed
        .combine(EXTEND_TARGETS.out.bed)
        // combine bed files into single list
        .map{ meta, bed1, bed2 ->
            [ meta, bed1 + [ bed2 ] ]
        }

    //
    // MODULE: Run bedtools/merge
    //
    def ch_bedtools_merge_extend_input = EXTEND_TARGETS.out.bed
        // add meta object to bed file
        .map{ bed ->
            [ [ id: bed.simpleName ], bed ]
        }
    BEDTOOLS_MERGE_EXTEND(ch_bedtools_merge_extend_input)
    ch_position_count_input = ch_position_count_input
        // add bedtools merge output bed
        .combine(BEDTOOLS_MERGE_EXTEND.out.bed
            // remove meta object from bed file
            .map{ meta, bed ->
                bed
            }
        )
        // combine bed files into single list
        .map{ meta, bed1, bed2 ->
            [ meta, bed1 + [ bed2 ] ]
        }

    //
    // MODULE: Run split_targets
    //
    SPLIT_TARGETS(BEDTOOLS_MERGE_EXTEND.out.bed)
    ch_position_count_input = ch_position_count_input
        // add split targets output bed
        .combine(SPLIT_TARGETS.out.bed
            // remove meta object from bed file
            .map{ meta, bed ->
                bed
            }
            )
        // combine bed files into single list
        .map{ meta, bed1, bed2 ->
            [ meta, bed1 + [ bed2 ] ]
        }

    //
    // MODULE: Run position_count
    //
    POSITION_COUNT(ch_position_count_input, FAI_TO_CHR_SIZES.out.tsv)

    //
    // MODULE: Run gatk4/preprocessintervals
    //
    GATK4_PREPROCESSINTERVALS(ch_reference, ch_reference_fai, ch_reference_dict, SPLIT_TARGETS.out.bed, channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/addorreplacereadgroups
    //
    GATK4_ADDORREPLACEREADGROUPS(ch_samplesheet, channel.value([ [], [] ]), channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/collectreadcounts
    //
    def ch_gatk4_collectreadcounts_input = GATK4_ADDORREPLACEREADGROUPS.out.bam
        // join gatk4 addorreplacereadgroups output bam and bai according to meta object
        .join(GATK4_ADDORREPLACEREADGROUPS.out.bai)
        // add gatk4 preprocessintervals output interval list
        .combine(
            GATK4_PREPROCESSINTERVALS.out.interval_list
                // remove meta object from interval list
                .map{ meta, list ->
                    list
                }
        )
    GATK4_COLLECTREADCOUNTS(ch_gatk4_collectreadcounts_input, channel.value([ [], [] ]), channel.value([ [], [] ]), channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/createreadcountpanelofnormals
    //
    def ch_gatk4_createreadcountpanelofnormals_input = GATK4_COLLECTREADCOUNTS.out.hdf5
        // remove meta object from hdf5 file
        .map{ meta, hdf5 ->
            hdf5
        }
        // combine all hdf5 files in channel into single list
        .collect()
        // add meta object to list of hdf5 files, id is file count + '_sample_pon'
        .map{ it ->
            [ [ id: it.size() + '_sample_pon' ], it ]
        }
    GATK4_CREATEREADCOUNTPANELOFNORMALS(ch_gatk4_createreadcountpanelofnormals_input)

    //
    // MODULE: Run gatk4/denoisereadcounts
    //
    GATK4_DENOISEREADCOUNTS(GATK4_COLLECTREADCOUNTS.out.hdf5, GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon)

    //
    // MODULE: Run gatk4/haplotypecaller
    //
    def ch_gatk4_haplotypecaller_input = ch_gatk4_collectreadcounts_input
        // add empty list
        .map{ meta, bam, bai, list ->
            [ meta, bam, bai, list, [] ]
        }
    GATK4_HAPLOTYPECALLER(ch_gatk4_haplotypecaller_input, ch_reference,  ch_reference_fai, ch_reference_dict, channel.value([ [], [] ]), channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/genomicsdbimport
    //
    def ch_gatk4_genomicsdbimport_input = GATK4_HAPLOTYPECALLER.out.vcf
        // create a tab-delimited file mapping of meta id to VCF file
        .map{ meta, vcf ->
            "${meta.id}\t${vcf.getName()}"
        }
        // write all mappings to a single file
        .collectFile(
            name: 'vcf_sample_map.txt',
            newLine: true
        )
        // concatenate gatk4 haplotypecaller output vcf
        .concat(
            GATK4_HAPLOTYPECALLER.out.vcf
                // remove meta object from vcf file
                .map{ meta, vcf ->
                    vcf
                }
        )
        // combine mapping file and vcf files in channel into single list
        .collect()
        // wrap list into another list
        .map{ it ->
            [ it ]
        }
        // add gatk4 haplotypecaller output tbi
        .combine(
            GATK4_HAPLOTYPECALLER.out.tbi
                // remove meta object from tbi file
                .map{ meta, tbi ->
                    tbi
                }
                // combine all tbi files in channel into single list
                .collect()
                // wrap list into another list
                .map{ tbi ->
                    [ tbi ]
                }
        )
        // add split targets output bed
        .combine(
            SPLIT_TARGETS.out.bed
                // remove meta object from bed file
                .map{ meta, bed ->
                    bed
                }
        )
        // add meta object and two empty lists to list of vcf/tbi/bed files, id is file count + '_sample_pon'
        .map{ vcf, tbi, bed ->
            [ [ id: tbi.size() + '_sample_pon' ], vcf, tbi, bed, [], [] ]
        }
    GATK4_GENOMICSDBIMPORT(ch_gatk4_genomicsdbimport_input, channel.value(false), channel.value(false), channel.value(true))

    //
    // MODULE: Run gatk4/genotypegvcfs
    //
    def ch_gatk4_gentypegvcfs_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb
        // add three empty lists to list of genomicsdb files
        .map{ meta, db ->
            [ meta, db, [], [], [] ]
        }
    GATK4_GENOTYPEGVCFS(ch_gatk4_gentypegvcfs_input, ch_reference, ch_reference_fai, ch_reference_dict, channel.value([ [], [] ]), channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/variantrecalibrator
    //
    def ch_gatk4_variantrecalibrator_input_1 = GATK4_GENOTYPEGVCFS.out.vcf
        // add gatk4 genotypegvcfs output tbi according to meta object
        .join(GATK4_GENOTYPEGVCFS.out.tbi)
    def ch_gatk4_variantrecalibrator_input_2 = channel.fromPath(gatk4_resources + '/*.vcf{,.gz}')
        // combine all vcf files in channel into single list
        .collect()
    def ch_gatk4_variantrecalibrator_input_3 = channel.fromPath(gatk4_resources + '/*.tbi')
        // combine all tbi files in channel into single list
        .collect()
    GATK4_VARIANTRECALIBRATOR(ch_gatk4_variantrecalibrator_input_1, ch_gatk4_variantrecalibrator_input_2, ch_gatk4_variantrecalibrator_input_3, channel.value(gatk4_resource_labels.split(';')), ch_reference_simple, ch_reference_fai_simple, ch_reference_dict_simple)

    //
    // MODULE: Run gatk4/applyvqsr
    //
    def ch_gatk4_applyvqsr_input = ch_gatk4_variantrecalibrator_input_1
        // join gatk4 variantrecalibrator output recal according to meta object
        .join(GATK4_VARIANTRECALIBRATOR.out.recal)
        // join gatk4 variantrecalibrator output idx according to meta object
        .join(GATK4_VARIANTRECALIBRATOR.out.idx)
        // join gatk4 variantrecalibrator output tranches according to meta object
        .join(GATK4_VARIANTRECALIBRATOR.out.tranches)
    GATK4_APPLYVQSR(ch_gatk4_applyvqsr_input, ch_reference_simple, ch_reference_fai_simple, ch_reference_dict_simple)

    //
    // MODULE: Run bcftools/norm
    //
    def ch_bcftools_norm_input = GATK4_APPLYVQSR.out.vcf
        // join gatk4 applyvqsr output tbi according to meta object
        .join(GATK4_APPLYVQSR.out.tbi)
    BCFTOOLS_NORM(ch_bcftools_norm_input, ch_reference)

    //
    // MODULE: Run extract_filtered_gatk_calls
    //
    EXTRACT_FILTERED_GATK_CALLS(BCFTOOLS_NORM.out.vcf)

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
