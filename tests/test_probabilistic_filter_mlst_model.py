"""Tests for the ProbabilisticFilterMlstSchemeModel class."""

# pylint: disable=redefined-outer-name

from pathlib import Path
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from xspect.handlers.pubmlst import PubMLSTHandler
from xspect.models.probabilistic_filter_mlst_model import (
    ProbabilisticFilterMlstSchemeModel,
)
from xspect.definitions import get_xspect_model_path

handler = PubMLSTHandler()
ORGANISM = "abaumannii"
SCHEME = "MLST (Oxford)"
SCHEME_ID = "1"


@pytest.fixture
def filter_model(tmpdir):
    """Fixture for the ProbabilisticFilterMlstSchemeModel class."""
    base_path = tmpdir.mkdir("xspect_data")
    return ProbabilisticFilterMlstSchemeModel(
        k=21,
        model_display_name=SCHEME,
        base_path=Path(base_path),
        scheme_url=f"https://rest.pubmlst.org/db/pubmlst_{ORGANISM}_seqdef/schemes/{SCHEME_ID}",
        organism=ORGANISM,
    )


@pytest.fixture
def trained_filter_model(filter_model, tmpdir):
    """Fixture for the ProbabilisticFilterModel class with trained model."""
    allele_path = Path(tmpdir) / "alleles"
    allele_path.mkdir()

    handler.download_alleles(ORGANISM, SCHEME, allele_path)

    filter_model.fit(allele_path)
    return filter_model


def test_model_initialization(
    filter_model,
):
    """Test the initialization of the ProbabilisticFilterMlstSchemeModel class."""
    assert filter_model.k == 21
    assert filter_model.model_display_name == SCHEME
    assert filter_model.organism == ORGANISM
    assert filter_model.model_type == "MLST"
    assert filter_model.fpr == 0.001


def test_model_save(trained_filter_model):
    """Test the save method of the ProbabilisticFilterMlstSchemeModel class."""
    trained_filter_model.save()
    assert (
        trained_filter_model.base_path / "abaumannii-mlst-oxford-mlst.json"
    ).exists()


def test_fit(trained_filter_model):
    """Test the fit method of the ProbabilisticFilterMlstSchemeModel class."""
    assert len(trained_filter_model.indices) == 7  # Amount of cobs_structures
    expected_values = {
        "Oxf_cpn60": 265,
        "Oxf_gdhB": 381,
        "Oxf_gltA": 241,
        "Oxf_gpi": 541,
        "Oxf_gyrB": 352,
        "Oxf_recA": 254,
        "Oxf_rpoD": 286,
    }
    for locus, size in expected_values.items():
        # Important: size can be greater, because the database is updated regularly
        assert trained_filter_model.loci.get(locus) >= size


def test_predict(trained_filter_model):
    """Test the predict method of the ProbabilisticFilterMlstSchemeModel class."""
    # Allele_ID_4 of Oxf_cpn60 with 401 kmers of length 21 each
    seq = Seq(
        "ATGAACCCAATGGATTTAAAACGCGGTATCGACATTGCAGTAAAAACTGTAGTTGAAAAT"
        "ATCCGTTCTATTGCTAAACCAGCTGATGATTTCAAAGCAATTGAACAAGTAGGTTCAATC"
        "TCTGCTAACTCTGATACTACTGTTGGTAAACTTATTGCTCAAGCAATGGAAAAAGTAGGT"
        "AAAGAAGGCGTAATCACTGTAGAAGAAGGTTCTGGCTTCGAAGACGCATTAGACGTTGTA"
        "GAAGGTATGCAGTTTGACCGTGGTTATATCTCTCCGTACTTTGCAAACAAACAAGATACT"
        "TTAACTGCTGAACTTGAAAATCCGTTCATTCTTCTTGTTGATAAAAAAATCAGCAACATT"
        "CGTGAATTGATTTCTGTTTTAGAAGCAGTTGCTAAAACTGGTAAACCACTTCTTATCATC"
        "G"
    )
    seq_record = SeqRecord(seq)
    res = trained_filter_model.predict(seq_record).hits.get("test")[0]
    allele_id = res.get("Strain type").get("Oxf_cpn60")

    assert allele_id.get("Allele_ID_4") == 401


def test_model_load(trained_filter_model):
    """Test the load method of the ProbabilisticFilterMlstSchemeModel class."""
    trained_filter_model.save()
    loaded_model = ProbabilisticFilterMlstSchemeModel.load(
        trained_filter_model.base_path / "abaumannii-mlst-oxford-mlst.json"
    )
    assert loaded_model.k == 21
    assert loaded_model.model_display_name == SCHEME
    assert loaded_model.organism == ORGANISM
    assert loaded_model.model_type == "MLST"
    assert len(loaded_model.indices) == 7
    expected_values = {
        "Oxf_cpn60": 265,
        "Oxf_gdhB": 381,
        "Oxf_gltA": 241,
        "Oxf_gpi": 541,
        "Oxf_gyrB": 352,
        "Oxf_recA": 254,
        "Oxf_rpoD": 286,
    }
    for locus, size in expected_values.items():
        # Important: size can be greater, because the database is updated regularly
        assert loaded_model.loci.get(locus) >= size


def test_sequence_splitter():
    """Test the sequence split method on an arbitrary short sequence."""
    model = ProbabilisticFilterMlstSchemeModel(
        k=4,
        model_display_name="Test Filter",
        organism="Test Organism",
        base_path=get_xspect_model_path(),
        scheme_url="",
    )
    # len(seq) = 80; len(substring) = 20
    # k = 4 means each substring (except the first one) starts 3 (k - 1) base pairs earlier
    seq = "AGCTATTTCGCTGATGTCGACTGATCAAAAAGCCGGCGCGCTTTCGTATAGGCTAGCTACGACATACGATCGATCACTGA"
    res = model.sequence_splitter(seq, 20)
    # Does not assert to 4 because of 3 additional base pairs when sliced
    # (Last slice has 12 base pairs)
    assert len(res) == 5


def test_has_sufficient_score(trained_filter_model):
    """Tests if the kmer scores of a scheme are sufficient."""
    sufficient_dict = {
        "Scores": {
            "Oxf_cpn60": 265,
            "Oxf_gdhB": 381,
            "Oxf_gltA": 241,
            "Oxf_gpi": 541,
            "Oxf_gyrB": 352,
            "Oxf_recA": 254,
            "Oxf_rpoD": 286,
        }
    }

    not_sufficient_dict = {
        "Scores": {
            "Oxf_cpn60": 100,
            "Oxf_gdhB": 20,
            "Oxf_gltA": 6,
            "Oxf_gpi": 55,
            "Oxf_gyrB": 21,
            "Oxf_recA": 0,
            "Oxf_rpoD": 7,
        }
    }

    assert trained_filter_model.has_sufficient_score(
        sufficient_dict, trained_filter_model.avg_locus_bp_size
    )

    assert not trained_filter_model.has_sufficient_score(
        not_sufficient_dict, trained_filter_model.avg_locus_bp_size
    )
