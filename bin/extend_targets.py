"""Expanding input BED regions based on supplied sample coverage profiles."""
import argparse
import logging
import math

# Set up logging. The logging level is set to INFO, and the log messages will include the timestamp, log level, and message.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Daniel Vodak
version_string = "1.0.0:26-08-14"

# To do?
# - Make sure that region starts and ends don't get outside the chromosome bounds
#   (not happening with our modest extensions of exonic targets, but could be an issue for more general use).

arg_parser = argparse.ArgumentParser(
  description="Expand input BED regions in both directions. "
              "The expansion is done one base at a time, as long as >= 'min_depth' coverage is maintained "
              "among the specified fraction of input samples "
              "(these are represented by their coverage profiles). "
              "(Specified male and female sample counts are utilized instead of total sample counts "
              "for coverage considerations during expansion of regions on Y and X chromosomes.)")

arg_parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s " + version_string)

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

arg_dict = vars(arg_parser.parse_args())

input_targets = arg_dict["input_targets_bed"]
output_targets = arg_dict["output_targets_bed"]
coverage_tsv_list = arg_dict["coverage_tsv"]
min_depth = arg_dict["min_depth"]
max_extension = arg_dict["max_extension"]
sample_fraction = arg_dict["sample_fraction"]
male_sample_count = arg_dict["male_sample_count"]
female_sample_count = arg_dict["female_sample_count"]

total_sample_count = len(coverage_tsv_list)

logger.info("Number of input coverage files/number of considered samples: {}".format(total_sample_count))

# Check selected parameter values. In case of unreasonable values, exit with a helpful message.
if ((male_sample_count + female_sample_count) != total_sample_count):
  logger.error("The total sample count ({}) does not equal the sum of specified"
               " male sample count ({}) and female sample count ({}) values. Exiting.".format(total_sample_count, male_sample_count, female_sample_count))
  exit(1)

if ((sample_fraction < 0.1) or (sample_fraction > 1.0)):
  logger.error("Please specify a \"sample_fraction\" between 0.1 and 1. Exiting.")
  exit(2)

min_covered_auto = math.ceil(sample_fraction * total_sample_count)
min_covered_chrX = math.ceil(sample_fraction * female_sample_count)
min_covered_chrY = math.ceil(sample_fraction * male_sample_count)

logger.info("Minimum required number of covering samples (number of samples with read depth >= {} at given position):"
            " {} samples for autosome regions; {} samples for chrX regions; {} samples for chrY regions.".format(min_depth, min_covered_auto, min_covered_chrX, min_covered_chrY))

original_targets = []
position_to_coverages = {}

# load the input target data
with open(input_targets, "r") as original_regions_file:
    for line in original_regions_file:
        line_s = line.strip().split("\t")

        # convert the start and end coordinates to integers
        region_entry = [line_s[0], int(line_s[1]), int(line_s[2]), line_s[3]]

        original_targets.append(region_entry)

logger.info("Number of input target regions: " + str(len(original_targets)))

# load the input coverage data
counter = 0
for coverage_tsv in coverage_tsv_list:
    counter += 1
    logger.info("Processing coverage file #" + str(counter) + "/" + str(len(coverage_tsv_list)) + "..")
    with open(coverage_tsv, "r") as coverage_file:
        for line in coverage_file:
            line_s = line.strip().split("\t")

            chr = line_s[0]
            pos = int(line_s[1]) + int(line_s[3]) - 1
            cov = int(line_s[4])

            if chr not in position_to_coverages:
                position_to_coverages[chr] = {}
            if pos not in position_to_coverages[chr]:
                position_to_coverages[chr][pos] = []
            position_to_coverages[chr][pos].append(cov)

# expand each input target based on the supplied coverage data
# write the expanded targets into the output file
with open(output_targets, "w") as grown_file:
    counter = 0
    for original_target in original_targets:
        counter += 1
        logger.info("Processing target #" + str(counter) + "/" + str(len(original_targets)) + "..")

        # parse the original target information
        chrom = original_target[0]
        start = original_target[1]
        end = original_target[2]
        name = original_target[3]

        # skip targets located on chromosomes not mentioned in the input coverage files - issue a warning for each target
        if (chrom not in position_to_coverages):
            logger.warning("Target chromosome is not reflected in the supplied coverage files. The following target will be skipped: ".format(":".join(original_target)))
            continue

        # initiate the start and end values for the extended target
        new_start = start
        new_end = end

        # determine how far the start position can be shifted without the coverage dropping below the specified threshold
        # note: start positions are considered to be inside the BED region
        for shift in range(1, max_extension + 1):
            test_start = start - shift
            if test_start in position_to_coverages[chrom]:
                coverages = position_to_coverages[chrom][test_start]
                covered_samples = len([coval for coval in coverages if coval >= min_depth])

                if ((chrom == "chrX") and (covered_samples >= min_covered_chrX)):
                    new_start = test_start
                elif ((chrom == "chrY") and (covered_samples >= min_covered_chrY)):
                    new_start = test_start
                elif ((chrom not in ["chrX", "chrY"]) and (covered_samples >= min_covered_auto)):
                    new_start = test_start
                # stop extending in case of coverage drop
                else:
                    break
            # stop if there is no data available for the tested position
            else:
                break

        # determine how far the end position can be shifted without the coverage dropping below the specified threshold
        # note: end positions are considered to be outside the BED region
        for shift in range(0, max_extension):
            test_end = end + shift
            if test_end in position_to_coverages[chrom]:
                coverages = position_to_coverages[chrom][test_end]
                covered_samples = len([coval for coval in coverages if coval >= min_depth])

                if ((chrom == "chrX") and (covered_samples >= min_covered_chrX)):
                    new_end = test_end + 1
                elif ((chrom == "chrY") and (covered_samples >= min_covered_chrY)):
                    new_end = test_end + 1
                elif ((chrom not in ["chrX", "chrY"]) and (covered_samples >= min_covered_auto)):
                    new_end = test_end + 1
                # stop extending in case of coverage drop
                else:
                    break
            # stop if there is no data available for the tested position
            else:
                break

        # output the adjusted target
        grown_file.write("\t".join([chrom, str(new_start), str(new_end), name + "_grown"]) + "\n")