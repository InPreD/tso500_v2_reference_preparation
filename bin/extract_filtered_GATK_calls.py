#!/usr/bin/env python3
"""Parsing input multi-sample VCF, extracting only PASS SNV/MNV (germline) calls with sufficient support in >=n samples. 
Saving the output in both VCF and BED-like TSV formats."""

import argparse
import logging
import gzip
import sys

# Set up logging. The logging level is set to INFO, and the log messages will include the timestamp, log level, and message.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():
    arg_parser = argparse.ArgumentParser(
        description="Filter the input VCF file, save the selected variants in output VCF and BED-like TSV files.\n"
        "The filtering:\n"
        "- Variants with missing 'AS_FilterStatus' or 'AC' info tags are discarded.\n"
        "- Variants with INFO tag 'AS_FilterStatus' value other than PASS are discarded.\n"
        "- Variants with AD tag missing in the format field are discarded.\n"
        "- Variants not reaching <min_read_count> supporting reads in >= <min_sample_count> samples are discarded.\n"
        "- Variants on contigs other than chr1-chr22, chrX and chrY are discarded.\n"
        "- Variants other than SNVs and MNVs of length 2 are discarded.\n"
        "VCF output: The original header and format details are preserved. The following changes are applied:\n"
        "- The sample-wise 'AC' format field values are changed to './.' if there are fewer than <min_read_count> supporting reads for given variant in given sample.\n"
        "BED-like TSV output: Contains the same variants as the VCF output, but only selected details. The file has the following columns:\n"
        "(1) chromosome, (2) 0-based start position, (3) one-based end position, (4) variant ID in format 'chrom:start_pos:ref>alt', "
        "(5) ref. allele, (6) alt. allele, (7) allelic count observed in the dataset population"
    )

    # [required arguments]
    arg_parser.add_argument("-i", "--input_vcf", required=True,
                            help="Path to the input VCF file.")
    arg_parser.add_argument("-o", "--output_tsv", required=True,
                            help="Path to the output TSV file.")
    arg_parser.add_argument("-p", "--output_vcf", required=True,
                            help="Path to the output VCF file.")
    # [optional arguments]
    arg_parser.add_argument("-r", "--min_read_count", default=5,
                            help="The number of variant-supporting reads required for considering a variant present in given sample (values >0 expected).")
    arg_parser.add_argument("-s", "--min_sample_count", default=5,
                            help="Only PASS variants present in at least this many samples will be extracted (values >0 expected).")

    arg_dict = vars(arg_parser.parse_args())

    min_read_count = arg_dict["min_read_count"]
    min_sample_count = arg_dict["min_sample_count"]
    input_vcf_path = arg_dict["input_vcf"]
    output_tsv_path = arg_dict["output_tsv"]
    output_vcf_path = arg_dict["output_vcf"]

    # Check selected parameter values. In case of unreasonable values, exit with a helpful message.
    if min_read_count < 1 or min_sample_count < 1:
        logger.error("Please make sure that supplied 'min_read_count' and 'min_sample_count' values are > 0. Exiting.")
        sys.exit(1)

    with gzip.open(input_vcf_path, "rt") as ivcf_file, \
        open(output_tsv_path, "w") as otsv_file, \
        open(output_vcf_path, "w") as ovcf_file:

        for line in ivcf_file:
            line_s = line.strip().split("\t")

            # print the header lines to the VCF output without change
            if (line_s[0][0] == "#"):
                ovcf_file.write(line)
            else:
                chrom = line_s[0]

                # skip variants outside the autosomes and sex chromosomes
                if (chrom not in ["chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9",  "chr10",
                                  "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                                  "chr21", "chr22", "chrX",  "chrY"]):
                    continue
                
                pos_one = line_s[1]
                pos_zero = str(int(line_s[1]) - 1)
                ref = line_s[3]
                alts = line_s[4].split(",")
                info_fields = line_s[7].split(";")
                format_tags = line_s[8].split(":")

                # skip variants with missing AD format-field tag
                if "AD" not in format_tags:
                    continue

                AD_index = format_tags.index("AD")

                # create a dictionary with INFO field items that have 'key=value' structure
                info_dict = {item.split("=")[0]:item.split("=")[1] for item in info_fields if "=" in item}

                # skip variants with missing 'AC' and/or 'AS_FilterStatus' INFO field items
                if "AS_FilterStatus" not in info_dict or "AC" not in info_dict:
                    continue

                AS_FilterStatus_values = info_dict["AS_FilterStatus"].split(",")
                allelic_counts = info_dict["AC"].split(",")

                # processing individual alt alleles
                for alt_index, alt in enumerate(alts):
                    filter_value = AS_FilterStatus_values[alt_index]
                    allelic_count = allelic_counts[alt_index]

                    # skip variants with 'AS_FilterStatus' value other than PASS
                    # keep only SNVs and MNVs with length <= 2
                    if filter_value != "PASS" or len(ref) != len(alt) or len(ref) > 2:
                        continue

                    # build the output VCF line, starting with information that is not sample-specific (columns 1-9)
                    ovcf_line = line_s[:9]
                    sample_count = 0
                    # check the sample-wise values, append potentially altered information
                    for sample_field in line_s[9:]:
                        sample_format_values = sample_field.split(":")
                        AD_output = sample_format_values[AD_index]
                        alt_AD_value = AD_output.split(",")[alt_index + 1]
                        if (int(alt_AD_value) < min_read_count):
                            AD_output = "./."
                        else:
                            sample_count += 1
                        ovcf_line.append(":".join(sample_format_values[:AD_index] + [AD_output] + sample_format_values[AD_index+1:]))

                    # output information for variants present in enough samples
                    if (sample_count >= min_sample_count):
                        ovcf_file.write("\t".join(ovcf_line) + "\n")

                        var_id = f"{chrom}:{pos_one}:{ref}:{alt}"
                        otsv_file.write("\t".join([chrom, pos_zero, pos_one, var_id, ref, alt, allelic_count]) + "\n")


if __name__ == "__main__":
    main()