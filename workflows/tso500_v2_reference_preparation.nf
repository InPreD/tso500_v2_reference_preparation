/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BEDTOOLS_BAMTOBED                         } from '../modules/nf-core/bedtools/bamtobed'
include { BEDTOOLS_COVERAGE                         } from '../modules/nf-core/bedtools/coverage'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_SUBTRACT } from '../modules/nf-core/bedtools/merge'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_EXTEND   } from '../modules/nf-core/bedtools/merge'
include { BEDTOOLS_SLOP                             } from '../modules/nf-core/bedtools/slop'
include { BEDTOOLS_SUBTRACT                         } from '../modules/nf-core/bedtools/subtract'
include { EXTEND_TARGETS                            } from '../modules/local/extend_targets'
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

    def ch_versions = channel.empty()
    .view()
    //
    // MODULE: Run fai to chromosome sizes
    //
    ch_fai_to_chr_sizes_input = channel.fromPath(reference_fai).map{ fai -> [ [ id: fai.simpleName ], fai ] }.first()
    FAI_TO_CHR_SIZES(ch_fai_to_chr_sizes_input)

    //
    // MODULE: Run sort
    //
    ch_sort_input = channel.fromPath(panel_bed).map{ bed -> [ [ id: bed.simpleName ], bed ] }
    SORT_BED(ch_sort_input)

    //
    // MODULE: Run bedtools/slop
    //
    BEDTOOLS_SLOP(SORT_BED.out.bed, FAI_TO_CHR_SIZES.out.tsv)
    ch_position_count_input = SORT_BED.out.bed.combine(BEDTOOLS_SLOP.out.bed.map{ meta, bed -> bed }).map{ meta, bed1, bed2 -> [ meta, [ bed1, bed2 ] ] }

    //
    // MODULE: Run bedtools/subtract
    //
    ch_bedtools_subtract_input = BEDTOOLS_SLOP.out.bed.map{ meta, bed -> [ meta, bed, file(panel_bed) ] }
    BEDTOOLS_SUBTRACT(ch_bedtools_subtract_input)
    ch_position_count_input = ch_position_count_input.combine(BEDTOOLS_SUBTRACT.out.bed.map{ meta, bed -> bed }).map{ meta, bed1, bed2 -> [ meta, bed1 + [ bed2 ] ] }

    //
    // MODULE: Run sort
    //
    SORT_SUBTRACT(BEDTOOLS_SUBTRACT.out.bed)

    //
    // MODULE: Run bedtools/merge
    //
    BEDTOOLS_MERGE_SUBTRACT(SORT_SUBTRACT.out.bed)
    ch_position_count_input = ch_position_count_input.combine(BEDTOOLS_MERGE_SUBTRACT.out.bed.map{ meta, bed -> bed }).map{ meta, bed1, bed2 -> [ meta, bed1 + [ bed2 ] ] }

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
    ch_bedtools_coverage_input = SORT_BAMTOBED.out.bed.combine(BEDTOOLS_MERGE_SUBTRACT.out.bed).map{ meta1, bed1, meta2, bed2 -> [ [ id: meta1.id + '_' + meta2.id, sex: meta1.sex.toLowerCase() ], bed2, bed1 ] }
    BEDTOOLS_COVERAGE(ch_bedtools_coverage_input, FAI_TO_CHR_SIZES.out.sorted_tsv.first())

    //
    // MODULE: Run extend_targets
    //
    ch_extend_targets_input = BEDTOOLS_COVERAGE.out.bed.map{it -> [ it ]}.collect().map{
        // count sex occurences and combine all bed files
        it ->
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
    ch_position_count_input = ch_position_count_input.combine(EXTEND_TARGETS.out.bed).map{ meta, bed1, bed2 -> [ meta, bed1 + [ bed2 ] ] }

    //
    // MODULE: Run bedtools/merge
    //
    ch_bedtools_merge_extend_input = EXTEND_TARGETS.out.bed.map{ bed -> [ [ id: bed.simpleName ], bed ] }
    BEDTOOLS_MERGE_EXTEND(ch_bedtools_merge_extend_input)
    ch_position_count_input = ch_position_count_input.combine(BEDTOOLS_MERGE_EXTEND.out.bed.map{ meta, bed -> bed }).map{ meta, bed1, bed2 -> [ meta, bed1 + [ bed2 ] ] }

    //
    // MODULE: Run split_targets
    //
    SPLIT_TARGETS(BEDTOOLS_MERGE_EXTEND.out.bed)
    ch_position_count_input = ch_position_count_input.combine(SPLIT_TARGETS.out.bed.map{ meta, bed -> bed }).map{ meta, bed1, bed2 -> [ meta, bed1 + [ bed2 ] ] }

    //
    // MODULE: Run split_targets
    //
    POSITION_COUNT(ch_position_count_input, FAI_TO_CHR_SIZES.out.tsv)

    //
    // MODULE: Run gatk4/preprocessintervals
    //
    ch_gatk4_preprocessintervals_input_1 = channel.fromPath(reference).map{ fa -> [ [ id: fa.simpleName ], fa ] }.first()
    ch_gatk4_preprocessintervals_input_2 = channel.fromPath(reference_dict).map{ dict -> [ [ id: dict.simpleName ], dict ] }.first()
    GATK4_PREPROCESSINTERVALS(ch_gatk4_preprocessintervals_input_1, ch_fai_to_chr_sizes_input, ch_gatk4_preprocessintervals_input_2, SPLIT_TARGETS.out.bed, channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/addorreplacereadgroups
    //
    GATK4_ADDORREPLACEREADGROUPS(ch_samplesheet, channel.value([ [], [] ]), channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/collectreadcounts
    //
    ch_gatk4_collectreadcounts_input = GATK4_ADDORREPLACEREADGROUPS.out.bam.join(GATK4_ADDORREPLACEREADGROUPS.out.bai).combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map{ meta, list -> list })
    GATK4_COLLECTREADCOUNTS(ch_gatk4_collectreadcounts_input, channel.value([ [], [] ]), channel.value([ [], [] ]), channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/createreadcountpanelofnormals
    //
    ch_gatk4_createreadcountpanelofnormals_input = GATK4_COLLECTREADCOUNTS.out.hdf5.map{ meta, hdf5 -> hdf5 }.collect().map{ it -> [ [ id: it.size() + '_sample_pon' ], it ]}
    GATK4_CREATEREADCOUNTPANELOFNORMALS(ch_gatk4_createreadcountpanelofnormals_input)

    //
    // MODULE: Run gatk4/denoisereadcounts
    //
    GATK4_DENOISEREADCOUNTS(GATK4_COLLECTREADCOUNTS.out.hdf5, GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon)

    //
    // MODULE: Run gatk4/haplotypecaller
    //
    ch_gatk4_haplotypecaller_input = ch_gatk4_collectreadcounts_input.map{ meta, bam, bai, list -> [ meta, bam, bai, list, [] ]}
    GATK4_HAPLOTYPECALLER(ch_gatk4_haplotypecaller_input, ch_gatk4_preprocessintervals_input_1,  ch_fai_to_chr_sizes_input, ch_gatk4_preprocessintervals_input_2, channel.value([ [], [] ]), channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/genomicsdbimport
    //
    ch_gatk4_genomicsdbimport_input = GATK4_HAPLOTYPECALLER.out.vcf.map{
        meta, vcf -> "${meta.id}\t${vcf.getName()}"
    }.collectFile(name: 'vcf_sample_map.txt', newLine: true).concat(GATK4_HAPLOTYPECALLER.out.vcf.map{meta, vcf -> vcf}).collect().map{
        it -> [ it ] 
    }.combine(GATK4_HAPLOTYPECALLER.out.tbi.map{ meta, tbi -> tbi }.collect().map{ tbi -> [ tbi ]}).combine(SPLIT_TARGETS.out.bed.map{meta, bed -> bed}).map{
        vcf, tbi, bed ->
        [ [ id: tbi.size() + '_sample_pon' ], vcf, tbi, bed, [], [] ]
    }
    GATK4_GENOMICSDBIMPORT(ch_gatk4_genomicsdbimport_input, channel.value(false), channel.value(false), channel.value(true))

    //
    // MODULE: Run gatk4/genotypegvcfs
    //
    ch_gatk4_gentypegvcfs_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, db -> [ meta, db, [], [], [] ] }
    GATK4_GENOTYPEGVCFS(ch_gatk4_gentypegvcfs_input, ch_gatk4_preprocessintervals_input_1, ch_fai_to_chr_sizes_input, ch_gatk4_preprocessintervals_input_2, channel.value([ [], [] ]), channel.value([ [], [] ]))

    //
    // MODULE: Run gatk4/variantrecalibrator
    //
    ch_gatk4_variantrecalibrator_input_1 = GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi)
    ch_gatk4_variantrecalibrator_input_2 = channel.fromPath(gatk4_resources + '/*.vcf???').collect()
    ch_gatk4_variantrecalibrator_input_3 = channel.fromPath(gatk4_resources + '/*.tbi').collect()
    GATK4_VARIANTRECALIBRATOR(ch_gatk4_variantrecalibrator_input_1, ch_gatk4_variantrecalibrator_input_2, ch_gatk4_variantrecalibrator_input_3, channel.value(gatk4_resource_labels.split(';')), ch_gatk4_preprocessintervals_input_1.map{ meta, fasta -> fasta }, ch_fai_to_chr_sizes_input.map{ meta, fai -> fai }, ch_gatk4_preprocessintervals_input_2.map{ meta, dict -> dict })

    //
    // MODULE: Run gatk4/applyvqsr
    //
    ch_gatk4_applyvqsr_input = ch_gatk4_variantrecalibrator_input_1.join(GATK4_VARIANTRECALIBRATOR.out.recal).join(GATK4_VARIANTRECALIBRATOR.out.idx).join(GATK4_VARIANTRECALIBRATOR.out.tranches)
    GATK4_APPLYVQSR(ch_gatk4_applyvqsr_input, ch_gatk4_preprocessintervals_input_1.map{ meta, fasta -> fasta }, ch_fai_to_chr_sizes_input.map{ meta, fai -> fai }, ch_gatk4_preprocessintervals_input_2.map{ meta, dict -> dict })

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
