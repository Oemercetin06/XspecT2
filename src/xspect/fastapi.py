"""FastAPI application for XspecT."""

import datetime
from pathlib import Path
from shutil import copyfileobj
from fastapi import FastAPI, UploadFile, BackgroundTasks
from xspect.definitions import get_xspect_runs_path, get_xspect_upload_path
from xspect.download_filters import download_test_filters
import xspect.model_management as mm
from xspect.models.result import StepType
from xspect.pipeline import ModelExecution, Pipeline, PipelineStep
from xspect.train import train_ncbi

app = FastAPI()


@app.get("/download-filters")
def download_filters():
    """Download filters."""
    download_test_filters("https://xspect2.s3.eu-central-1.amazonaws.com/models.zip")


@app.get("/classify")
def classify(genus: str, file: str, meta: bool = False, step: int = 500):
    """Classify uploaded sample."""

    path = get_xspect_upload_path() / file

    pipeline = Pipeline(genus + " classification", "Test Author", "test@example.com")
    species_execution = ModelExecution(
        genus.lower() + "-species", sparse_sampling_step=step
    )
    if meta:
        species_filtering_step = PipelineStep(
            StepType.FILTERING, genus, 0.7, species_execution
        )
        genus_execution = ModelExecution(
            genus.lower() + "-genus", sparse_sampling_step=step
        )
        genus_execution.add_pipeline_step(species_filtering_step)
        pipeline.add_pipeline_step(genus_execution)
    else:
        pipeline.add_pipeline_step(species_execution)

    run = pipeline.run(Path(path))
    time_str = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    save_path = get_xspect_runs_path() / f"run_{time_str}.json"
    run.save(save_path)

    return run.to_dict()


@app.post("/train")
def train(genus: str, background_tasks: BackgroundTasks, svm_steps: int = 1):
    """Train NCBI model."""
    background_tasks.add_task(train_ncbi, genus, svm_steps)

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
