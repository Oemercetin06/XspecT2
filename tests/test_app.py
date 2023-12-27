"""Test the WebApp."""

# pylint: disable=redefined-outer-name

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
        assert session["quick"] == True
        assert session["metagenome"] == False
        assert session["OXA"] == False
        assert "GCF_000240185.1_ASM24018v2_genomic.fna.txt" in session["filename"]
        assert "success" in response.text


def test_post_ic(client):
    """Test the clonetype assignment (ClassT) session config."""
    with client:
        response = client.post(
            "/ic",
            json=[
                "TTGATCGGTGCGTTGGCAACAAAAAAAT",
                "GCF_000240185.1_ASM24018v2_genomic.fna",
                True,  # quick
                True,  # IC 1
                True,  # IC 2
                True,  # IC 3
                True,  # IC 4
                True,  # IC 5
                True,  # IC 6
                True,  # IC 7
                True,  # IC 8
                False,  # Added
                False,  # OXA
            ],
        )
        assert session["quick"] == True
        # this includes mystery "added" checkbox
        assert session["IC_lookup"] == [True] * 8 + [False]
        assert session["OXA"] == False
        assert "GCF_000240185.1_ASM24018v2_genomic.fna.txt" in session["filename"]
        assert "success" in response.text


def test_species_results(client):
    """Test the species (Xspect) assignment & result page."""
    # get file
    with client.session_transaction() as session:
        session[
            "filename"
        ] = "files/f2f522083ec7c705GCF_000240185.1_ASM24018v2_genomic.fna.txt"
        session["quick"] = True
        session["metagenome"] = False
        session["OXA"] = False
        session["genus"] = "Acinetobacter"

    response = client.get("/assignspec", follow_redirects=True)
    assert len(response.history) == 1
    assert response.request.path == "/resultsspec"


def test_ic_results(client):
    """Test the clonetype assignment (ClassT) result page."""
    # get file
    with client.session_transaction() as session:
        session[
            "filename"
        ] = "files/f2f522083ec7c705GCF_000240185.1_ASM24018v2_genomic.fna.txt"
        session["quick"] = True
        session["IC_lookup"] = [True] * 8 + [False]
        session["OXA"] = False

    response = client.get("/assign", follow_redirects=True)
    assert len(response.history) == 1
    assert response.request.path == "/results"
