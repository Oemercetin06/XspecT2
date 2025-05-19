"""FastAPI-based web application for XspecT."""

from uuid import uuid4
import json
from shutil import copyfileobj
import importlib.resources as pkg_resources
from fastapi import APIRouter, FastAPI, HTTPException, UploadFile, BackgroundTasks
from fastapi.responses import RedirectResponse
from xspect.definitions import get_xspect_runs_path, get_xspect_upload_path
from xspect.download_models import download_test_models
import xspect.model_management as mm
from xspect.train import train_from_ncbi
from xspect import classify, filter_sequences
from fastapi.staticfiles import StaticFiles

app = FastAPI()
app.mount(
    "/xspect-web",
    StaticFiles(directory=str(pkg_resources.files("xspect") / "xspect-web" / "dist")),
    name="static",
)
router = APIRouter()


@app.get("/")
def root():
    """Root endpoint, forwards to /xspect-web/index.html."""
    return RedirectResponse(url="/xspect-web/index.html")


@router.get("/download-filters")
def download_filters():
    """Download filters."""
    download_test_models("http://assets.adrianromberg.com/xspect-models.zip")


@router.get("/classification-result")
def get_classification_result(uuid: str):
    """Get classification result."""
    result_path = get_xspect_runs_path() / f"result_{uuid}.json"
    if not result_path.exists():
        raise HTTPException(
            status_code=404, detail="No result found for the specified uuid."
        )
    return json.loads(result_path.read_text())


@router.post("/classify")
def classify_post(
    classification_type: str,
    model: str,
    file: str,
    background_tasks: BackgroundTasks,
    step: int = 1,
):
    """Classify uploaded sample."""
    input_path = get_xspect_upload_path() / file
    if not input_path.exists():
        raise FileNotFoundError(f"File {input_path} does not exist.")

    uuid = str(uuid4())

    if classification_type == "Genus":
        background_tasks.add_task(
            classify.classify_genus,
            model,
            input_path,
            get_xspect_runs_path() / f"result_{uuid}.json",
            step=step,
        )
        return {"message": "Classification started.", "uuid": uuid}

    elif classification_type == "Species":
        background_tasks.add_task(
            classify.classify_species,
            model,
            input_path,
            get_xspect_runs_path() / f"result_{uuid}.json",
            step=step,
        )
        return {"message": "Classification started.", "uuid": uuid}

    raise NotImplementedError(
        f"Classification type {classification_type} is not implemented."
    )


router.post("/filter")


def filter_post(
    filter_type: str,
    model: str,
    input_file: str,
    threshold: float,
    filter_species: str = None,
):
    """Filter sequences."""
    input_path = get_xspect_upload_path() / input_file
    output_path = get_xspect_upload_path() / f"filtered_{input_file}"

    if not input_path.exists():
        raise FileNotFoundError(f"File {input_path} does not exist.")

    if filter_type == "Genus":
        filter_sequences.filter_genus(model, input_path, output_path, threshold)
        return {"message": "Genus Filtering started."}

    elif filter_type == "Species":
        filter_sequences.filter_species(
            model, filter_species, input_path, output_path, threshold
        )
        return {"message": "Species Filtering started."}

    raise NotImplementedError(f"Filter type {filter_type} is not implemented.")


@router.post("/train")
def train(genus: str, background_tasks: BackgroundTasks, svm_steps: int = 1):
    """Train NCBI model."""
    background_tasks.add_task(train_from_ncbi, genus, svm_steps)

    return {"message": "Training started."}


@router.get("/list-models")
def list_models():
    """List available models."""
    return mm.get_models()


@router.get("/model-metadata")
def get_model_metadata(model_slug: str):
    """Get metadata of a model."""
    return mm.get_model_metadata(model_slug)


@router.post("/model-metadata")
def post_model_metadata(model_slug: str, author: str, author_email: str):
    """Update metadata of a model."""
    try:
        mm.update_model_metadata(model_slug, author, author_email)
    except ValueError as e:
        return {"error": str(e)}
    return {"message": "Metadata updated."}


@router.post("/model-display-name")
def post_model_display_name(model_slug: str, filter_id: str, display_name: str):
    """Update display name of a filter in a model."""
    try:
        mm.update_model_display_name(model_slug, filter_id, display_name)
    except ValueError as e:
        return {"error": str(e)}
    return {"message": "Display name updated."}


@router.post("/upload-file")
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


app.include_router(router, prefix="/api", tags=["api"])
