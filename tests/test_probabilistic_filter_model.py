"""Tests for the ProbabilisticFilterModel class."""

# pylint: disable=redefined-outer-name, line-too-long, protected-access

from pathlib import Path
import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from xspect.models.probabilistic_filter_model import ProbabilisticFilterModel


@pytest.fixture
def filter_model(tmpdir):
    """Fixture for the ProbabilisticFilterModel class."""
    base_path = tmpdir.mkdir("xspect_data")
    return ProbabilisticFilterModel(
        k=21,
        model_display_name="Test Filter",
        author="John Doe",
        author_email="john.doe@example.com",
        model_type="Species",
        base_path=Path(base_path),
    )


@pytest.fixture
def trained_filter_model(filter_model, multiple_assembly_dir_path):
    """Fixture for the ProbabilisticFilterModel class with trained model."""
    filter_model.fit(Path(multiple_assembly_dir_path))
    return filter_model


def test_model_initialization(filter_model):
    """Test the initialization of the ProbabilisticFilterModel class."""
    assert filter_model.k == 21
    assert filter_model.model_display_name == "Test Filter"
    assert filter_model.author == "John Doe"
    assert filter_model.author_email == "john.doe@example.com"
    assert filter_model.model_type == "Species"


def test_fit(filter_model, multiple_assembly_dir_path):
    """Test the fit method of the ProbabilisticFilterModel class."""
    filter_model.fit(Path(multiple_assembly_dir_path))
    cobs_path = filter_model.base_path / filter_model.slug() / "index.cobs_classic"
    assert cobs_path.exists()
    assert filter_model.display_names == {
        "GCF_000006945": "GCF_000006945.2_ASM694v2_genomic",
        "GCF_000018445": "GCF_000018445.1_ASM1844v1_genomic",
        "GCF_000069245": "GCF_000069245.1_ASM6924v1_genomic",
    }


def test_fit_with_display_names(filter_model, multiple_assembly_dir_path):
    """Test the fit method of the ProbabilisticFilterModel class with display names."""
    display_names = {
        "GCF_000006945.2_ASM694v2_genomic": "Salmonella enterica",
        "GCF_000018445.1_ASM1844v1_genomic": "Acinetobacter baumannii ACICU",
        "GCF_000069245.1_ASM6924v1_genomic": "Acinetobacter baumannii AYE",
    }
    # cobs only uses the names until the first "." internally
    expected_display_names = {
        "GCF_000006945": "Salmonella enterica",
        "GCF_000018445": "Acinetobacter baumannii ACICU",
        "GCF_000069245": "Acinetobacter baumannii AYE",
    }
    filter_model.fit(Path(multiple_assembly_dir_path), display_names=display_names)
    cobs_path = filter_model.base_path / filter_model.slug() / "index.cobs_classic"
    assert cobs_path.exists()
    assert filter_model.display_names == expected_display_names


def test_predict(trained_filter_model):
    """Test the predict method of the ProbabilisticFilterModel class."""

    salmonella_record = SeqRecord(
        seq=Seq(
            "AGAGATTACGTCTGGTTGCAAGAGATCATGACAGGGGGAATTGGTTGAAAATAAATATATCGCCAGCAGCACATGAACAA"
        ),
        id="test",
    )

    res = trained_filter_model.predict(salmonella_record)
    assert res.get_scores()["total"] == {
        "GCF_000006945": 1.0,
        "GCF_000069245": 0.0,
        "GCF_000018445": 0.0,
    }
    assert res.get_total_hits() == {
        "GCF_000006945": 60,
        "GCF_000069245": 0,
        "GCF_000018445": 0,
    }


@pytest.mark.parametrize(
    ["assembly_file_path", "expected_scores"],
    [
        (
            "GCF_000069245.1_ASM6924v1_genomic.fna",
            {"GCF_000006945": 0.01, "GCF_000018445": 0.63, "GCF_000069245": 1.0},
        ),
        (
            "GCF_000018445.1_ASM1844v1_genomic.fna",
            {"GCF_000006945": 0.01, "GCF_000018445": 1.0, "GCF_000069245": 0.64},
        ),
        (
            "GCF_000006945.2_ASM694v2_genomic.fna",
            {"GCF_000006945": 1.0, "GCF_000018445": 0.0, "GCF_000069245": 0.0},
        ),
    ],
    indirect=["assembly_file_path"],
)
def test_predict_from_file(trained_filter_model, assembly_file_path, expected_scores):
    """Test the predict method of the ProbabilisticFilterModel class with a file
    containing multiple (sub)sequences."""
    res = trained_filter_model.predict(Path(assembly_file_path))
    assert res.get_scores()["total"] == expected_scores


def test_model_save_and_load(trained_filter_model):
    """Test the save and load methods of the ProbabilisticFilterModel class."""
    trained_filter_model.save()
    loaded_model = ProbabilisticFilterModel.load(
        trained_filter_model.base_path / "test-filter-species.json"
    )
    assert trained_filter_model.to_dict() == loaded_model.to_dict()


def test_trained_filter_model_save_and_load(trained_filter_model):
    """Test the save and load methods of the ProbabilisticFilterModel class."""
    trained_filter_model.save()
    loaded_model = ProbabilisticFilterModel.load(
        trained_filter_model.base_path / "test-filter-species.json"
    )
    assert trained_filter_model.to_dict() == loaded_model.to_dict()


def test_count_kmers(filter_model):
    """Test the _count_kmers method of the ProbabilisticFilterModel class."""
    sequence = Seq(
        "AGAGATTACGTCTGGTTGCAAGAGATCATGACAGGGGGAATTGGTTGAAAATAAATATATCGCCAGCAGCACATGAACAA"
    )
    kmer_counts = filter_model._count_kmers(sequence)
    # k = 21
    assert kmer_counts == 60


def test_calculate_hits(trained_filter_model):
    """Test the calculate_hits method of the ProbabilisticFilterModel class."""
    salmonella_sequence = Seq(
        "AGAGATTACGTCTGGTTGCAAGAGATCATGACAGGGGGAATTGGTTGAAAATAAATATATCGCCAGCAGCACATGAACAA"
    )

    for step in range(1, 5):
        hits = trained_filter_model.calculate_hits(salmonella_sequence, step=step)
        assert hits == {
            "GCF_000006945": 60 / step,
            "GCF_000069245": 0,
            "GCF_000018445": 0,
        }
