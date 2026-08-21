"""Splitting input BED regions into fragments of specified size."""
import argparse
import logging
import math
from typing import List

# Set up logging. The logging level is set to INFO, and the log messages will include the timestamp, log level, and message.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# split given region into one or more fragments according to supplied size restrictions
def fragment_region(chrom: str, start: int, end: int, minimum_size: int, standard_size: int) -> List[List[str]]:
    fragments = []

    target_size = end - start

    # determine the number of fragments created during input target splitting
    # 1) non-overlapping fragments of selected standard size (as many as possible, n >= 0) are created (broken off from the input target), starting at the 'start' position
    # 2-a) if the input target remainder has fewer than 'minimum_size' bp in length (a "short remainder")
    #   2-a-a) and there is at least one fragment of standard size: the remainder region will be appended to the last fragment of standard size
    #   2-a-b) and there are no fragments of standard size: the short remainder will be output as the only fragment (it will be identical to the whole original target)
    # 2-b) if the target region remainer has more than 'minimum_size' bp in length (a "long remainder"), there will be an extra fragment created (its size will be between 'minimum_size' bp and 'standard_size' bp)
    # - the extra fragment will be identical to the whole original target if there were no standard size fragments
    standard_fragments = math.floor(target_size / standard_size)
    remainder_size = target_size % standard_size
    extra_fragment = False
    if (remainder_size >= minimum_size):
        extra_fragment = True

    # output the individual fragments
    # 1) - standard fragment output
    for fragment_index in range(standard_fragments):
        fragment_start = start + standard_size*fragment_index
        fragment_end = start + standard_size*(fragment_index + 1)

        # the 2-a-a) case: a "short remainder" extending the last standard fragment (fragment size > standard_size)
        if ((fragment_index == (standard_fragments - 1)) and (not extra_fragment)):
            fragment_end = end

        fragments.append([chrom, f"{fragment_start}", f"{fragment_end}"])

    # the 2-a-b) case: a "short remainder" with no standard fragments (input target size < minimum_size)
    if ((not extra_fragment) and (standard_fragments == 0)):
        fragments.append([chrom, f"{start}", f"{end}"])

    # the 2-b) case: a long remainder (minimum_size <= remainder length < standard_size) that will form an extra fragment
    if extra_fragment:
        fragment_start = start + standard_size*standard_fragments
        fragment_end = end

        fragments.append([chrom, f"{fragment_start}", f"{fragment_end}"])

    return(fragments)

def main():
    arg_parser = argparse.ArgumentParser(
    description="Split the input BED regions into output fragments of specified lengh ('standard_size' bp)."
                " The splitting won't allow creation of fragments shorter than 'minimum_size' bp"
                " (input regions already shorter than 'minimum_size' will be included in the output as-is)."
                " Depending on the setting, an output fragment can be up to 'standard_size'*2-1 bp long."
                " Please note: The output BED file will only contain columns 1-3. Each input region is handled independetly of the others."
                " Example 1: Input target of 540 bp would be split into intervals of 150 bp, 150 bp and 240 bp with the default settings."
                " Example 2: Input target of 560 bp would be split into intervals of 150 bp, 150 bp, 150 bp and 110 bp with the default settings.")

    arg_parser.add_argument("--input_targets_bed", required=True,
                            help="path to the original target BED file")
    arg_parser.add_argument("--output_targets_bed", required=True,
                            help="path to the expanded target BED file")
    arg_parser.add_argument("--minimum_size", type=int, default=100,
                            help="minimum size (in bp) of the output intervals")
    arg_parser.add_argument("--standard_size", type=int, default=150,
                            help="standard size (in bp) of the output intervals")
    arg_parser.add_argument("--output_name_column", type=bool, default=False,
                            help="whether the name column should be included in the output")

    arg_dict = vars(arg_parser.parse_args())

    input_targets = arg_dict["input_targets_bed"]
    output_targets = arg_dict["output_targets_bed"]
    minimum_size = arg_dict["minimum_size"]
    standard_size = arg_dict["standard_size"]
    output_name_column = arg_dict["output_name_column"]

    logger.info("BED region splitting initiated.")

    input_region_count = 0
    output_region_count = 0

    # Check selected parameter values. In case of unreasonable values, exit with a helpful message.
    if (minimum_size > standard_size):
        logger.error("Please make sure that standard_size is larger than minimum size. Exiting.")
        exit(1)

    if ((minimum_size < 1) or (standard_size < 1)):
        logger.error("Please make sure that the minimum_size and standard_size are larger than 0. Exiting.")
        exit(2)

    # Iterate through the input file one target region at a time, perform the necessary splitting
    # and output the newly created fragments into the output file
    with open(input_targets, "r") as input_file, \
        open(output_targets, "w") as output_file:

        for line in input_file:
            line_s = line.strip().split("\t")

            # skip incomplete input BED regions, issue a warning for each
            line_element_count = len(line_s)
            if (line_element_count < 3):
                logger.warning(f"Too few columns/fields on the following input BED file line (\"{line_s}\"). The line will be skipped.")
                continue

            input_region_count += 1
            fragment_count = 0

            # determine input target properties
            chrom = line_s[0]
            start = int(line_s[1])
            end = int(line_s[2])
            
            # determine target name; construct a new name if there is none available on input
            target_name = f"region_{input_region_count}"
            if (line_element_count > 3):
                target_name = line_s[3]

            # split the input region into fragments
            region_fragments = fragment_region(chrom, start, end, minimum_size, standard_size)

            # print the individual fragments into the output file
            for region_fragment in region_fragments:
                output_region_count += 1
                fragment_count += 1

                output_line = region_fragment
                if output_name_column:
                    output_line += [f"{target_name}|{fragment_count}"]

                output_file.write("\t".join(output_line) + "\n")

    logger.info("BED region splitting done.")
    logger.info(f"Number of processed input regions: {input_region_count}")
    logger.info(f"Number of output regions: {output_region_count}")

if __name__ == "__main__":
    main()