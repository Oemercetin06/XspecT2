"""
Tests for the FastAPI module.
"""

# pylint: disable=redefined-outer-name

import json
import time
import pytest
from fastapi.testclient import TestClient
from xspect.definitions import get_xspect_runs_path
from xspect.web import app
from xspect.model_management import get_model_metadata
from pathlib import Path
from Bio import SeqIO


@pytest.fixture
def client():
    """Create a FastAPI test client."""
    return TestClient(app)


@pytest.fixture
def client_with_uploaded_file(client, request):
    """Create a FastAPI test client with an uploaded file."""
    assembly_file_path = Path(request.param)
    with open(assembly_file_path, "rb") as f:
        client.post(
            "/api/upload-file",
            files={"file": (assembly_file_path.name, f)},
        )
    return client


def test_list_models(client):
    """Test the /list-models endpoint."""
    response = client.get("/api/list-models")
    assert response.status_code == 200
    assert "Genus" in response.json()
    assert "Species" in response.json()


def test_get_model_metadata(client):
    """Test the /model-metadata endpoint."""
    response = client.get(
        "/api/model-metadata", params={"model_slug": "acinetobacter-species"}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["model_display_name"] == "Acinetobacter"
    assert response_json["display_names"]["471"] == "Acinetobacter calcoaceticus"


def test_post_model_metadata(client):
    """Test the /model-metadata endpoint."""
    response = client.post(
        "/api/model-metadata",
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
        "/api/model-display-name",
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
    response = client_with_uploaded_file.post(
        "/api/classify",
        params={
            "classification_type": "Species",
            "model": "Acinetobacter",
            "file": Path(assembly_file_path).name,
            "step": 1,
        },
    )
    assert response.status_code == 200
    uuid = response.json()["uuid"]
    result_file_path = Path(get_xspect_runs_path()) / f"result_{uuid}.json"

    # wait for background task to finish. Check every 5 seconds for 1 minute
    for _ in range(12):
        if result_file_path.exists():
            break
        time.sleep(5)

    assert result_file_path.exists()
    result = json.loads(result_file_path.read_text())
    assert result["prediction"] == "470"


@pytest.mark.parametrize(
    ["assembly_file_path", "client_with_uploaded_file", "filter_type", "species"],
    [
        (
            "tests/test_assemblies/GCF_000018445.1_ASM1844v1_genomic.fna",
            "tests/test_assemblies/GCF_000018445.1_ASM1844v1_genomic.fna",
            "Genus",
            None,
        ),
        (
            "tests/test_assemblies/GCF_000018445.1_ASM1844v1_genomic.fna",
            "tests/test_assemblies/GCF_000018445.1_ASM1844v1_genomic.fna",
            "Species",
            "470",
        ),
    ],
    indirect=["client_with_uploaded_file"],
)
def test_filter(client_with_uploaded_file, assembly_file_path, filter_type, species):
    """Test the /filter endpoint."""
    response = client_with_uploaded_file.post(
        "/api/filter",
        params={
            "filter_type": filter_type,
            "genus": "Acinetobacter",
            "input_file": Path(assembly_file_path).name,
            "threshold": 0.7,
            "filter_species": species,
        },
    )
    assert response.status_code == 200
    uuid = response.json()["uuid"]

    result_file_path = Path(get_xspect_runs_path()) / f"result_{uuid}.json"
    filtered_file_path = Path(get_xspect_runs_path()) / f"filtered_{uuid}.fasta"

    # wait for background task to finish. Check every 5 seconds for 1 minute
    for _ in range(12):
        if result_file_path.exists() and filtered_file_path.exists():
            break
        time.sleep(5)

    assert result_file_path.exists()
    assert filtered_file_path.exists()

    response = client_with_uploaded_file.get(f"/api/filtering-result?uuid={uuid}")
    assert response.status_code == 200
    assert json.loads(response.text)["uuid"] == uuid

    input_sequences = SeqIO.parse(assembly_file_path, "fasta")
    output_sequences = SeqIO.parse(filtered_file_path, "fasta")

    # the chosen threshold means that the output should contain fewer sequences than the input
    assert len(list(input_sequences)) > len(list(output_sequences))

    response = client_with_uploaded_file.get(f"/api/download-filtered?uuid={uuid}")
    assert response.status_code == 200
    assert (
        response.headers["Content-Disposition"]
        == f'attachment; filename="filtered_{uuid}.fasta"'
    )
