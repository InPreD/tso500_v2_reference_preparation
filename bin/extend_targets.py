"""Expanding input BED regions based on supplied sample coverage profiles."""
import argparse
import logging
import math
from typing import Dict, List

# extend given region based on
# 1) the supplied coverage profiles (here reduced to the number of samples with 'high-enough' coverage at the relevant positions)
# 2) the required number of samples with 'high-enough' coverage in the extension areas 
# 3) the selected extension limit
def extend_region(chrom: str, start: int, end: int, relevant_position_dict: Dict[int, int], max_extension: int, coverage_threshold: int) -> List[str]:

    # initiate the start and end values for the extended region
    new_start = start
    new_end = end

    # determine how far the start position can be shifted without an unacceptable coverage drop across the coverage profiles
    # note: start positions are considered to be inside the BED region
    for shift in range(1, max_extension + 1):
        test_start = start - shift
        if test_start in relevant_position_dict:
            covered_samples = relevant_position_dict[test_start]

            if (covered_samples >= coverage_threshold):
                new_start = test_start
            # stop extending in case of a coverage drop
            else:
                break
        # stop if there is no data available for the tested position
        else:
            break

    # determine how far the end position can be shifted without an unacceptable coverage drop across the coverage profiles
    # note: end positions are considered to be outside the BED region
    for shift in range(0, max_extension):
        test_end = end + shift
        if test_end in relevant_position_dict:
            covered_samples = relevant_position_dict[test_end]

            if (covered_samples >= coverage_threshold):
                new_end = test_end + 1
            # stop extending in case of a coverage drop
            else:
                break
        # stop if there is no data available for the tested position
        else:
            break

    # return the extended region as a List of str elements
    return([chrom, f"{new_start}", f"{new_end}"])

# Set up logging. The logging level is set to INFO, and the log messages will include the timestamp, log level, and message.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# To do:
# - Make sure that region starts and ends don't get outside the chromosome bounds
#   (not happening with our modest extensions of exonic targets, but could be an issue for more general use).
# - Allow supplying coverage files via a text file with paths.

arg_parser = argparse.ArgumentParser(
  description="Expand input BED regions in both directions. "
              "The expansion is done one base at a time, as long as >= 'min_depth' coverage is maintained "
              "among the specified fraction of input samples "
              "(these are represented by their coverage profiles). "
              "(Specified male and female sample counts are utilized instead of total sample counts "
              "for coverage considerations during expansion of regions on Y and X chromosomes.)")

arg_parser.add_argument("--input_targets_bed", required=True,
                        help="path to the input BED file")
arg_parser.add_argument("--output_targets_bed", required=True,
                        help="path to the output (expanded) BED file")
arg_parser.add_argument("--coverage_tsv", action="append", required=True,
                        help="path to a reference sample coverage TSV file with 5 columns, no header: "
                             "1. chromosome; 2. region-start (inclusive, 0-based); 3. region-end (non-inclusive, 0-based); "
                             "4. position index within region (starting from 1); 5. coverage value")
arg_parser.add_argument("--min_depth", type=int, default=50,
                        help="required minimum required among the reference samples")
arg_parser.add_argument("--max_extension", type=int, default=250,
                        help="maximum allowed target extension in each direction")
arg_parser.add_argument("--sample_fraction", type=float, default=0.7,
                        help="the fraction of samples that need to have the required minimum coverage")
arg_parser.add_argument("--male_sample_count", type=int, required=True,
                        help="number of male patient samples")
arg_parser.add_argument("--female_sample_count", type=int, required=True,
                        help="number of female patient samples")
arg_parser.add_argument("--output_name_column", type=bool, default=True,
                        help="whether the name column should be included in the output")

arg_dict = vars(arg_parser.parse_args())

input_targets = arg_dict["input_targets_bed"]
output_targets = arg_dict["output_targets_bed"]
coverage_tsv_list = arg_dict["coverage_tsv"]
min_depth = arg_dict["min_depth"]
max_extension = arg_dict["max_extension"]
sample_fraction = arg_dict["sample_fraction"]
male_sample_count = arg_dict["male_sample_count"]
female_sample_count = arg_dict["female_sample_count"]
output_name_column = arg_dict["output_name_column"]

total_sample_count = len(coverage_tsv_list)

logger.info("BED region extension initiated.")
logger.info(f"Number of input coverage files/number of considered samples: {total_sample_count}")

# Check selected parameter values. In case of unreasonable values, exit with a helpful message.
if ((male_sample_count + female_sample_count) != total_sample_count):
  logger.error(f"The total sample count ({total_sample_count}) does not equal the sum of specified"
               " male sample count ({male_sample_count}) and female sample count ({female_sample_count}) values. Exiting.")
  exit(1)

if ((sample_fraction < 0.1) or (sample_fraction > 1.0)):
  logger.error("Please specify a \"sample_fraction\" parameter value between 0.1 and 1. Exiting.")
  exit(2)

if ((min_depth < 10) or (max_extension < 10)):
  logger.error("Please specify \"min_depth\" and \"max_extension\" parameter values larger than 9. Exiting.")
  exit(3)

# determine how many samples with high-enough coverage (depth >= min_depth) at given position
# will be required for including that position in the extended regions (the output)
# - autosomes, chrX and chrY have all their own threshold values
min_covered_auto = math.ceil(sample_fraction * total_sample_count)
min_covered_chrX = math.ceil(sample_fraction * female_sample_count)
min_covered_chrY = math.ceil(sample_fraction * male_sample_count)

logger.info(f"Minimum required number of covering samples (number of samples with read depth >= {min_depth} at given position):"
            " {min_covered_auto} samples for autosome regions; {min_covered_chrX} samples for chrX regions; {min_covered_chrY} samples for chrY regions.")

# variables for holding the input data
original_targets = []
position_to_covered_samples = {}

# load the input target data
with open(input_targets, "r") as original_regions_file:
    for line in original_regions_file:
        line_s = line.strip().split("\t")

        # skip incomplete input BED regions, issue a warning for each
        line_element_count = len(line_s)
        if (line_element_count < 3):
            logger.warning(f"Too few columns/fields on the following input BED file line (\"{line_s}\"). The line will be skipped.")
            continue

        # determine target name; construct a new name if there is none available on input
        target_name = f"region_{len(original_targets) + 1}"
        if (line_element_count > 3):
            target_name = line_s[3]

        # construct a region record
        region_entry = [line_s[0], line_s[1], line_s[2], target_name]
        original_targets.append(region_entry)

logger.info(f"Number of input target regions: {len(original_targets)}")

# load the input coverage data
file_counter = 0
for coverage_tsv in coverage_tsv_list:
    file_counter += 1
    logger.info(f"Processing coverage file #{file_counter}/{len(coverage_tsv_list)}..")
    with open(coverage_tsv, "r") as coverage_file:
        for line in coverage_file:
            line_s = line.strip().split("\t")

            chr = line_s[0]
            pos = int(line_s[1]) + int(line_s[3]) - 1
            cov = int(line_s[4])

            if (cov >= min_depth):
                if chr not in position_to_covered_samples:
                    position_to_covered_samples[chr] = {}
                if pos not in position_to_covered_samples[chr]:
                    position_to_covered_samples[chr][pos] = 1
                else:
                    position_to_covered_samples[chr][pos] += 1

# expand each input target based on the supplied coverage data
# write the expanded targets into the output file
region_count = 0
with open(output_targets, "w") as output_file:
    for original_target in original_targets:
        region_count += 1
        logger.info(f"Processing target #{region_count}/{len(original_targets)}..")

        # parse the original target information
        chrom = original_target[0]
        start = int(original_target[1])
        end = int(original_target[2])
        name = original_target[3]

        # skip targets located on chromosomes not mentioned in the input coverage files - issue a warning for each target
        if (chrom not in position_to_covered_samples):
            logger.warning(f"The following target will be skipped, as its chromosome is not reflected in the supplied coverage files: {':'.join(original_target)}")
            continue

        # select the correct coverage threshold based on the chromosome
        coverage_threshold = min_covered_auto
        if (chrom == "chrX"):
            coverage_threshold = min_covered_chrX
        elif(chrom == "chrY"):
            coverage_threshold = min_covered_chrY

        # extract only the relevant parts of the 'position_to_covered_samples' dictionary
        # (only data for the positions that might be queried in the subsequent extension process for given region)
        relevant_position_list = [pos for pos in position_to_covered_samples[chrom] if (start - max_extension) <= pos <= (end + max_extension - 1)]
        relevant_position_dict = {pos: position_to_covered_samples[chrom][pos] for pos in relevant_position_list}

        # extend the input region
        extended_region = extend_region(chrom, start, end, relevant_position_dict, max_extension, coverage_threshold)

        output_line = extended_region
        if output_name_column:
            output_line += [f"{name}|extended"]

        output_file.write("\t".join(output_line) + "\n")

logger.info("BED region extension done.")
logger.info(f"Number of processed regions: {region_count}")