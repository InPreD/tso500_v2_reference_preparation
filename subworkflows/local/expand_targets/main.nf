//
// Subworkflow with functionality specific to the inpred/tso500_v2_reference_preparation pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_BAMTOBED                         } from '../../../modules/nf-core/bedtools/bamtobed'
include { BEDTOOLS_COVERAGE                         } from '../../../modules/nf-core/bedtools/coverage'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_SUBTRACT } from '../../../modules/nf-core/bedtools/merge'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_EXTEND   } from '../../../modules/nf-core/bedtools/merge'
include { BEDTOOLS_SLOP                             } from '../../../modules/nf-core/bedtools/slop'
include { BEDTOOLS_SUBTRACT                         } from '../../../modules/nf-core/bedtools/subtract'
include { EXTEND_TARGETS                            } from '../../../modules/local/extend_targets'
include { FAI_TO_CHR_SIZES                          } from '../../../modules/local/fai_to_chr_sizes'
include { POSITION_COUNT                            } from '../../../modules/local/position_count'
include { SORT as SORT_BED                          } from '../../../modules/local/sort'
include { SORT as SORT_SUBTRACT                     } from '../../../modules/local/sort'
include { SORT as SORT_BAMTOBED                     } from '../../../modules/local/sort'
include { SPLIT_TARGETS                             } from '../../../modules/local/split_targets'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO EXPAND TARGETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXPAND_TARGETS {

    take:
    ch_reference_fai // channel: reference fai file
    ch_samplesheet   // channel: samplesheet read in from --input
    panel_bed        // string: path to panel bed file

    main:

    // general channel
    def ch_versions = channel.empty()

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

    emit:
    bed         = SPLIT_TARGETS.out.bed
    versions    = ch_versions
}