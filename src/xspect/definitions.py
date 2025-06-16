"""This module contains definitions for the XspecT package."""

from pathlib import Path
from os import getcwd

fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]
fastq_endings = ["fastq", "fq"]


def get_xspect_root_path() -> Path:
    """
    Return the root path for XspecT data.

    Returns the path to the XspecT data directory, which can be located either in the user's home directory or in the current working directory.
    If neither exists, it creates the directory in the user's home directory.

    Returns:
        Path: The path to the XspecT data directory.
    """

    home_based_dir = Path.home() / "xspect-data"
    if home_based_dir.exists():
        return home_based_dir

    cwd_based_dir = Path(getcwd()) / "xspect-data"
    if cwd_based_dir.exists():
        return cwd_based_dir

    home_based_dir.mkdir(exist_ok=True, parents=True)
    return home_based_dir


def get_xspect_model_path() -> Path:
    """
    Return the path to the XspecT models.

    Returns the path to the XspecT models directory, which is located within the XspecT data directory.
    If the directory does not exist, it creates the directory.

    Returns:
        Path: The path to the XspecT models directory.
    """
    model_path = get_xspect_root_path() / "models"
    model_path.mkdir(exist_ok=True, parents=True)
    return model_path


def get_xspect_upload_path() -> Path:
    """
    Return the path to the XspecT upload directory.

    Returns the path to the XspecT uploads directory, which is located within the XspecT data directory.
    If the directory does not exist, it creates the directory.

    Returns:
        Path: The path to the XspecT uploads directory.
    """
    upload_path = get_xspect_root_path() / "uploads"
    upload_path.mkdir(exist_ok=True, parents=True)
    return upload_path


def get_xspect_runs_path() -> Path:
    """
    Return the path to the XspecT runs directory.

    Returns the path to the XspecT runs directory, which is located within the XspecT data directory.
    If the directory does not exist, it creates the directory.

    Returns:
        Path: The path to the XspecT runs directory.
    """
    runs_path = get_xspect_root_path() / "runs"
    runs_path.mkdir(exist_ok=True, parents=True)
    return runs_path


def get_xspect_mlst_path() -> Path:
    """
    Return the path to the XspecT MLST directory.

    Returns the path to the XspecT MLST directory, which is located within the XspecT data directory.
    If the directory does not exist, it creates the directory.

    Returns:
        Path: The path to the XspecT MLST directory.
    """
    mlst_path = get_xspect_root_path() / "mlst"
    mlst_path.mkdir(exist_ok=True, parents=True)
    return mlst_path
