from pytest import mark

from .split_targets import fragment_region

#fragment_region(chrom: str, start: int, end: int, minimum_size: int, standard_size: int) -> List[List[str]]

# the following parameter value limits are imposed by the python script that includes the tested function
# - minimum_size <= standard_size
# - minimum_size >= 1
# - standard_size >= 1
@mark.parametrize(
    "chrom, start, end, minimum_size, standard_size, fragment_list",
    [
        ("chr1", 100, 300, 75, 100, [["chr1", "100", "200"], ["chr1", "200", "300"]]), # all output fragments being of standard_size
        ("chr1", 100, 350, 75, 100, [["chr1", "100", "200"], ["chr1", "200", "350"]]), # the last fragment being larger than standard_size
        ("chr1", 100, 275, 50, 100, [["chr1", "100", "200"], ["chr1", "200", "275"]]), # the last fragment being smaller than standard_size (an "extra fragment")
        ("chr1", 100, 150, 75, 200, [["chr1", "100", "150"]]) # a fragment being smaller than minimum_size (the whole input region being preserved)
    ]
)

def test_fragment_region(chrom, start, end, minimum_size, standard_size, fragment_list):
    assert fragment_region(chrom, start, end, minimum_size, standard_size) == fragment_list
