from pathlib import Path
from xspect.handlers.pubmlst import PubMLSTHandler

handler = PubMLSTHandler()


def test_download_default(tmpdir):
    """Tests the download functionality of alleles for the Oxford MLST model (A. baumannii)."""
    allele_path = Path(tmpdir) / "oxford"
    handler.download_alleles("abaumannii", "MLST (Oxford)", allele_path)
    oxford_loci = [
        "Oxf_cpn60",
        "Oxf_gdhB",
        "Oxf_gltA",
        "Oxf_gpi",
        "Oxf_gyrB",
        "Oxf_recA",
        "Oxf_rpoD",
    ]
    for locus_path in sorted(allele_path.iterdir()):
        locus_name = locus_path.name
        assert locus_name in oxford_loci


def test_get_strain_type_name():
    """Tests the POST request that gets the strain type name based on the highest kmer results."""
    # url is from the oxford scheme of A.baumannii.
    post_url = "https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes/1"
    alleles = {  # Translates to ST 1 in the db.
        "Oxf_gyrB": 1,
        "Oxf_gltA": 1,
        "Oxf_gdhB": 1,
        "Oxf_recA": 1,
        "Oxf_gpi": 1,
        "Oxf_cpn60": 1,
        "Oxf_rpoD": 6,
    }

    strain_type_name = handler.get_strain_type_name(alleles, post_url)

    assert strain_type_name["ST"] == "1"
