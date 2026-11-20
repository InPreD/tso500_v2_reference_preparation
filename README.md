# inpred/tso500_v2_reference_preparation

## Introduction

**inpred/tso500_v2_reference_preparation** ingests normal-patient BAM files, the TSO500 panel BED file, and a reference genome with its index and sequence dictionary. It expands and splits panel targets, performs GATK preprocessing and variant recalibration, creates a panel of normals, and generates reference outputs for [tso500 vs post processing](https://github.com/InPreD/tso500_v2_post_processing/).

1. Expand panel targets and split intervals for analysis using bedtools and local Python utilities.
1. Collect read counts and allelic counts from normal BAM files using GATK4.
1. Generate a panel of normals and denoised read-count models using GATK4.
1. Filter heterozygous sites and extract reference intervals for downstream post-processing.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,sex,bam
sample_1,M,assets/test_stub/sample_1.bam
sample_2,F,assets/test_stub/sample_2.bam
```

Each row represents one normal sample:

- `sample`: Unique sample name.
- `sex`: Sample sex, using M or F.
- `bam`: Path to the sample BAM file.

Now, you can run the pipeline using:

```bash
nextflow run inpred/tso500_v2_reference_preparation \
    -profile docker \
    --input samplesheet.csv \
    --panel_bed panel.bed \
    --reference reference.fa \
    --reference_dict reference.fa.dict \
    --reference_fai reference.fa.fai \
    --gatk4_resources gatk4_resources \
    --gatk4_resource_labels \
      '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz;--resource:omni,known=false,training=true,truth=false,prior=12.0 omni.vcf.gz;--resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz;--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz' \
    --outdir results
```

For stubbing run:

```bash
nextflow run inpred/tso500_v2_reference_preparation \
    -profile test_stub,docker \
    --outdir test_stub_outdir
```

> [!IMPORTANT]
> Pipeline parameters can be provided via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

inpred/tso500_v2_reference_preparation was originally written by @marrip and @danielvo.

We thank the following people for their extensive assistance in the development of this pipeline:

- @Omer0191
- @Stautis

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).
