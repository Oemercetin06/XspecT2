import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from xspect.mlst_feature.mlst_helper import pick_scheme, pick_scheme_from_models_dir
from xspect.mlst_feature.pub_mlst_handler import PubMLSTHandler
from xspect.models.probabilistic_filter_mlst_model import (
    ProbabilisticFilterMlstSchemeModel,
)
from xspect.definitions import get_xspect_model_path

handler = PubMLSTHandler()
PubMLSTHandler.download_alleles(handler, False)


@pytest.fixture
def filter_model():
    """Fixture for the ProbabilisticFilterMlstSchemeModel class."""
    return ProbabilisticFilterMlstSchemeModel(
        k_value=21,
        model_name="Test Filter",
        base_path=get_xspect_model_path() / "test",
        scheme_url="https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes/1",
    )


@pytest.fixture
def trained_filter_model(filter_model, monkeypatch):
    """Fixture for the ProbabilisticFilterModel class with trained model."""
    monkeypatch.setattr("builtins.input", lambda _: "1")
    scheme = pick_scheme(handler.get_scheme_paths())
    filter_model.fit(scheme)
    return filter_model


def test_model_initialization(
    filter_model,
):
    """Test the initialization of the ProbabilisticFilterMlstSchemeModel class."""
    assert filter_model.k_value == 21
    assert filter_model.model_name == "Test Filter"
    assert filter_model.model_type == "Strain"
    assert filter_model.fpr == 0.001


def test_model_save(trained_filter_model):
    """Test the save method of the ProbabilisticFilterMlstSchemeModel class."""
    trained_filter_model.save()
    assert (
        trained_filter_model.base_path / "MLST (Oxford)" / "MLST (Oxford).json"
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


def test_predict(trained_filter_model, monkeypatch):
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
    monkeypatch.setattr("builtins.input", lambda _: "1")  # 1 = Oxford, 2 = Pasteur
    scheme = pick_scheme_from_models_dir()
    # [0] = Dict with the highest hits ([1] has all results which are not needed)
    res = trained_filter_model.predict(scheme, seq_record).hits.get("test")[0]
    allele_id = res.get("Strain type").get("Oxf_cpn60")

    assert allele_id.get("Allele_ID_4") == 401


def test_model_load(trained_filter_model):
    """Test the load method of the ProbabilisticFilterMlstSchemeModel class."""
    loaded_model = ProbabilisticFilterMlstSchemeModel.load(
        get_xspect_model_path() / "test" / "MLST" / "MLST (Oxford)"
    )
    assert loaded_model.k_value == 21
    assert loaded_model.model_name == "Test Filter"
    assert loaded_model.model_type == "Strain"
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
        k_value=4,
        model_name="Test Filter",
        base_path=get_xspect_model_path(),
        scheme_url="",
    )
    # len(seq) = 80; len(substring) = 20
    # k = 4 means each substring (except the first one) starts 3 (k - 1) base pairs earlier
    seq = "AGCTATTTCGCTGATGTCGACTGATCAAAAAGCCGGCGCGCTTTCGTATAGGCTAGCTACGACATACGATCGATCACTGA"
    res = model.sequence_splitter(seq, 20)
    # Does not assert 4 because of 3 additional base pairs when sliced (Last slice has 12 base pairs)
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
