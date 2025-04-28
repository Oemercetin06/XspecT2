"""FastAPI application for XspecT."""

from uuid import uuid4
from pathlib import Path
from shutil import copyfileobj
from fastapi import FastAPI, UploadFile, BackgroundTasks
from xspect.definitions import get_xspect_runs_path, get_xspect_upload_path
from xspect.download_models import download_test_models
from xspect.file_io import filter_sequences
import xspect.model_management as mm
from xspect.train import train_from_ncbi

app = FastAPI()


@app.get("/download-filters")
def download_filters():
    """Download filters."""
    download_test_models("http://assets.adrianromberg.com/xspect-models.zip")


@app.get("/classify")
def classify(genus: str, file: str, meta: bool = False, step: int = 500):
    """Classify uploaded sample."""

    input_path = get_xspect_upload_path() / file

    uuid = str(uuid4())

    if meta:
        genus_model = mm.get_genus_model(genus)
        genus_result = genus_model.predict(input_path, step=step)
        included_ids = genus_result.get_filtered_subsequence_labels(genus)
        if not included_ids:
            return {"message": "No sequences found for the given genus."}
        filtered_path = get_xspect_runs_path() / f"filtered_{uuid}.fasta"
        filter_sequences(
            Path(input_path),
            Path(filtered_path),
            included_ids=included_ids,
        )
        input_path = filtered_path

    species_model = mm.get_species_model(genus)
    species_result = species_model.predict(input_path, step=step)
    species_result.save(get_xspect_runs_path() / f"result_{uuid}.json")
    return species_result.to_dict()


@app.post("/train")
def train(genus: str, background_tasks: BackgroundTasks, svm_steps: int = 1):
    """Train NCBI model."""
    background_tasks.add_task(train_from_ncbi, genus, svm_steps)

    return {"message": "Training started."}


@app.get("/list-models")
def list_models():
    """List available models."""
    return mm.get_models()


@app.get("/model-metadata")
def get_model_metadata(model_slug: str):
    """Get metadata of a model."""
    return mm.get_model_metadata(model_slug)


@app.post("/model-metadata")
def post_model_metadata(model_slug: str, author: str, author_email: str):
    """Update metadata of a model."""
    try:
        mm.update_model_metadata(model_slug, author, author_email)
    except ValueError as e:
        return {"error": str(e)}
    return {"message": "Metadata updated."}


@app.post("/model-display-name")
def post_model_display_name(model_slug: str, filter_id: str, display_name: str):
    """Update display name of a filter in a model."""
    try:
        mm.update_model_display_name(model_slug, filter_id, display_name)
    except ValueError as e:
        return {"error": str(e)}
    return {"message": "Display name updated."}


@app.post("/upload-file")
def upload_file(file: UploadFile):
    """Upload file to the server."""
    upload_path = get_xspect_upload_path() / file.filename

    if not upload_path.exists():
        try:
            with upload_path.open("wb") as buffer:
                copyfileobj(file.file, buffer)
        finally:
            file.file.close()

    return {"filename": file.filename}
