"""This module contains definitions for the XspecT package."""

from pathlib import Path
from os import getcwd

fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]
fastq_endings = ["fastq", "fq"]


def get_xspect_root_path():
    """Return the root path for XspecT data."""
    root_path = Path(getcwd()) / "xspect-data"
    root_path.mkdir(exist_ok=True, parents=True)
    return root_path


def get_xspect_model_path():
    """Return the path to the XspecT models."""
    model_path = get_xspect_root_path() / "models"
    model_path.mkdir(exist_ok=True, parents=True)
    return model_path


def get_xspect_upload_path():
    """Return the path to the XspecT upload directory."""
    upload_path = get_xspect_root_path() / "uploads"
    upload_path.mkdir(exist_ok=True, parents=True)
    return upload_path


def get_xspect_runs_path():
    """Return the path to the XspecT runs directory."""
    runs_path = get_xspect_root_path() / "runs"
    runs_path.mkdir(exist_ok=True, parents=True)
    return runs_path


def get_xspect_mlst_path():
    """Return the path to the XspecT runs directory."""
    mlst_path = get_xspect_root_path() / "mlst"
    mlst_path.mkdir(exist_ok=True, parents=True)
    return mlst_path
