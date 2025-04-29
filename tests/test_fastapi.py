"""
Tests for the FastAPI module.
"""

# pylint: disable=redefined-outer-name

import pytest
from fastapi.testclient import TestClient
from xspect.fastapi import app
from xspect.model_management import get_model_metadata
from pathlib import Path


@pytest.fixture
def client():
    """Create a FastAPI test client."""
    return TestClient(app)


@pytest.fixture
def client_with_uploaded_file(client, request):
    """Create a FastAPI test client with an uploaded file."""
    assembly_file_path = Path(request.param)
    with open(assembly_file_path, "rb") as f:
        response = client.post(
            "/upload-file",
            files={"file": (assembly_file_path.name, f)},
        )
    return client


def test_list_models(client):
    """Test the /list-models endpoint."""
    response = client.get("/list-models")
    assert response.status_code == 200
    assert "Genus" in response.json()
    assert "Species" in response.json()


def test_get_model_metadata(client):
    """Test the /model-metadata endpoint."""
    response = client.get(
        "/model-metadata", params={"model_slug": "acinetobacter-species"}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["model_display_name"] == "Acinetobacter"
    assert response_json["display_names"]["471"] == "Acinetobacter calcoaceticus"


def test_post_model_metadata(client):
    """Test the /model-metadata endpoint."""
    response = client.post(
        "/model-metadata",
        params={
            "model_slug": "acinetobacter-species",
            "author": "Test Author",
            "author_email": "test@example.com",
        },
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["message"] == "Metadata updated."
    model_metadata = get_model_metadata("acinetobacter-species")
    assert model_metadata["author"] == "Test Author"
    assert model_metadata["author_email"] == "test@example.com"


def test_post_model_display_name(client):
    """Test the /model-display-name endpoint."""
    response = client.post(
        "/model-display-name",
        params={
            "model_slug": "acinetobacter-species",
            "filter_id": "470",
            "display_name": "AB",
        },
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["message"] == "Display name updated."
    model_metadata = get_model_metadata("acinetobacter-species")
    assert model_metadata["display_names"]["470"] == "AB"


@pytest.mark.parametrize(
    ["assembly_file_path", "client_with_uploaded_file"],
    [
        (
            "tests/test_assemblies/GCF_000018445.1_ASM1844v1_genomic.fna",
            "tests/test_assemblies/GCF_000018445.1_ASM1844v1_genomic.fna",
        )
    ],
    indirect=["client_with_uploaded_file"],
)
def test_classify(client_with_uploaded_file, assembly_file_path):
    """Test the /classify endpoint."""
    response = client_with_uploaded_file.get(
        "/classify",
        params={
            "genus": "acinetobacter",
            "file": Path(assembly_file_path).name,
            "step": 1,
            "included_ids": ["470"],
        },
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["prediction"] == "470"
