"""Test Flask web app."""

# pylint: disable=redefined-outer-name

import time
import pytest
from flask import session
import WebApp


@pytest.fixture()
def app():
    """Create test app."""
    yield WebApp.app


@pytest.fixture()
def client(app):
    """Create test client for the app."""
    return app.test_client()


@pytest.fixture(scope="module")
def genus_client(client, request):
    """Create test client with genus-specific filters already set up."""
    with client:
        client.post("/change_genus", data={"genus": request.param})
        return client


def test_home(client):
    """Test the home page."""
    response = client.get("/home")
    assert response.status_code == 200
    assert "Welcome to XspecT and ClAssT!" in response.text


def test_get_species(client):
    """Test the species (Xspect) page."""
    response = client.get("/species")
    assert response.status_code == 200
    assert (
        "Select Genus and upload Sequence Reads or a Genome Assembly" in response.text
    )


def test_get_ic(client):
    """Test the clonetype assignment (ClassT) page."""
    response = client.get("/ic")
    assert response.status_code == 200
    assert (
        "ClAssT is reasy to use! Upload data to identify International-Clones."
        in response.text
    )


def test_about(client):
    """Test the about page."""
    response = client.get("/about")
    assert response.status_code == 200
    assert "How-to-use" in response.text


def test_change_genus(client):
    """Test the change genus page."""
    with client:
        client.post("/change_genus", data={"genus": "Salmonella"})
        assert session["genus"] == "Salmonella"


def test_post_species(client):
    """Test the species (Xspect) session config."""
    with client:
        response = client.post(
            "/species",
            json=[
                "TTGATCGGTGCGTTGGCAACAAAAAAAT",
                "GCF_000240185.1_ASM24018v2_genomic.fna",
                True,  # quick
                False,  # metagenome
                False,  # OXA
            ],
        )
        assert session["quick"]
        assert not session["metagenome"]
        assert not session["OXA"]
        assert "GCF_000240185.1_ASM24018v2_genomic.fna.txt" in session["filename"]
        assert "success" in response.text


@pytest.mark.parametrize(
    ["genus_client", "assembly_file_path", "species", "oxa"],
    [
        (
            "Acinetobacter",
            "GCF_000069245.1_ASM6924v1_genomic.fna",
            "A 470",
            "[0.0, 0.0, 0.0, 0.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]",
        ),
        (
            "Acinetobacter",
            "GCF_000018445.1_ASM1844v1_genomic.fna",
            "A 470",
            "[0.0, 0.0, 0.0, 0.03, 0.0, 0.01, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0]",
        ),
        ("Salmonella", "GCF_000006945.2_ASM694v2_genomic.fna", "S 28901", ""),
    ],
    indirect=["genus_client", "assembly_file_path"],
)
def test_results(genus_client, assembly_file_path, species, oxa):
    """Test the species (Xspect) assignment & result page."""
    with genus_client.session_transaction() as session:
        session["filename"] = assembly_file_path
        session["quick"] = True
        session["metagenome"] = False
        session["OXA"] = bool(oxa)

    response = genus_client.get("/assignspec", follow_redirects=True)
    assert len(response.history) == 1
    assert response.request.path == "/resultsspec"
    assert species in response.text
    assert oxa in response.text
