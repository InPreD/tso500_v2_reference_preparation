from pytest import mark

from ..extend_targets import extend_region

#extend_region(chrom: str, start: int, end: int, relevant_position_dict: Dict[int, int], max_extension: int, coverage_threshold: int) -> List[str]

# the following parameter value limits are imposed by the python script that includes the tested function
# - max_extension >= 10
@mark.parametrize(
    "chrom, start, end, relevant_position_dict, max_extension, coverage_threshold, extended_region",
    [
        ("chr1", 10, 20, {9:12, 8:15, 7:10, 20:13, 21:15, 22:10}, 3, 10, ["chr1",  "7", "23"]), # full extension in both directions
        ("chr1", 10, 20, {9:12,       7:10, 20:13, 21:15, 22:10}, 3, 10, ["chr1",  "9", "23"]), # missing coverage data upstream, only partial upstream extension
        ("chr1", 10, 20, {9:12, 8:15, 7:10, 20:13,        22:10}, 3, 10, ["chr1",  "7", "21"]), # missing coverage data downstream, only partial downstream extension
        ("chr1", 10, 20, {9:9,  8:15, 7:10, 20:13, 21:15, 22:10}, 3, 10, ["chr1", "10", "23"]), # coverage < coverage_threshold upstream, no upstream extension
        ("chr1", 10, 20, {9:12, 8:15, 7:10, 20:8,  21:15, 22:8},  3, 10, ["chr1",  "7", "20"]), # coverage < coverage_threshold downstream, no downstream extension
        ("chr1", 10, 20, {      8:15, 7:10, 20:8,  21:15, 22:10}, 3, 10, ["chr1", "10", "20"]), # low coverage upstream & missing coverage data downstream, no extension
        ("chr1", 10, 20, {},                                      3, 10, ["chr1", "10", "20"])  # no coverage data, no extension
    ]
)

def test_extend_region(chrom, start, end, relevant_position_dict, max_extension, coverage_threshold, extended_region):
    assert extend_region(chrom, start, end, relevant_position_dict, max_extension, coverage_threshold) == extended_region
