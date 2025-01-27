from xspect.definitions import get_xspect_mlst_path
from xspect.mlst_feature.pub_mlst_handler import PubMLSTHandler

handler = PubMLSTHandler()


def test_download_default():
    """Tests the download functionality of alleles for the default case."""
    handler.download_alleles(False)
    species_path = get_xspect_mlst_path() / "abaumannii"
    first_scheme_path = species_path / "MLST (Oxford)"
    oxford_loci = [
        "Oxf_cpn60",
        "Oxf_gdhB",
        "Oxf_gltA",
        "Oxf_gpi",
        "Oxf_gyrB",
        "Oxf_recA",
        "Oxf_rpoD",
    ]
    assert sum(1 for _ in species_path.iterdir()) == 2  # 2 schemes
    for locus_path in sorted(first_scheme_path.iterdir()):
        locus_name = str(locus_path).split("/")[-1]
        assert locus_name in oxford_loci


def test_chosen_download(monkeypatch):
    """Tests the download functionality of alleles for a chosen scheme."""
    inputs = iter(["2", "1", "no"])  # Simulate input: [species, scheme, repeat]
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))
    handler.download_alleles(True)
    assert handler.scheme_list == [
        "https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes/1"
    ]
