"""Test the detection of misclassified sequences."""

from xspect.misclassification_detection.mapping import MappingHandler
from xspect.misclassification_detection.simulate_reads import extract_random_reads
from xspect.ncbi import NCBIHandler
import pytest, os

from xspect.misclassification_detection.point_pattern_analysis import (
    PointPatternAnalysis,
)


@pytest.fixture(scope="module")
def ncbi_handler():
    """Fixture for the NCBI class."""
    ncbi_api_key = os.environ.get("NCBI_API_KEY")
    return NCBIHandler(api_key=ncbi_api_key)


@pytest.fixture()
def reference(ncbi_handler, tmp_path):
    """Fixture for the path to the genome."""
    path = ncbi_handler.download_reference_genome(470, tmp_path)
    read_path = tmp_path / "reads.fasta"
    extract_random_reads(str(path), str(read_path), 150, 1000)
    return str(path), str(read_path)


@pytest.fixture()
def mapping_handler(reference):
    """Fixture for the path to the genome."""
    handler = MappingHandler(str(reference[0]), str(reference[1]))
    return handler


def test_extract_random_reads(reference, tmp_path):
    """Test the read simulation."""
    read_path = tmp_path / "read.fasta"
    extract_random_reads(str(reference[0]), str(read_path), 150, 10)
    assert read_path.exists()


def test_mapping_procedure(reference, mapping_handler):
    """Test the read simulation."""
    assert not os.path.isfile(mapping_handler.bam)
    mapping_handler.map_reads_onto_reference()
    mapping_handler.extract_starting_coordinates()
    genome_length = mapping_handler.get_total_genome_length()
    start_coordinates = mapping_handler.get_start_coordinates()
    assert genome_length > 3000000
    assert os.path.isfile(mapping_handler.bam)
    assert os.path.isfile(mapping_handler.tsv)
    assert len(start_coordinates) > 0


def test_point_pattern_analysis():
    """Test the point pattern density functions."""
    genome_length = 400000
    first_point_list = [  # Not clustered
        500,
        5000,
        10000,
        30000,
        80000,
        100000,
        104000,
        170000,
        200000,
        300000,
        399999,
    ]
    first_analysis = PointPatternAnalysis(first_point_list, genome_length)

    second_point_list = [  # clustered
        500,
        5000,
        10000,
        81000,
        81500,
        82000,
        84000,
        200000,
        300000,
        399999,
    ]
    second_analysis = PointPatternAnalysis(second_point_list, genome_length)

    assert first_analysis.ripleys_k()[0] == False
    assert first_analysis.ripleys_k_edge_corrected()[0] == False

    assert second_analysis.ripleys_k()[0] == True
    assert second_analysis.ripleys_k_edge_corrected()[0] == True
