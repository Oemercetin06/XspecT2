"""Download filters from public repository."""

import shutil
from tempfile import TemporaryDirectory
from pathlib import Path
import requests

from xspect.definitions import get_xspect_model_path


def download_test_models(url: str) -> None:
    """
    Download models from the specified URL.

    This function downloads a zip file from the given URL, extracts its contents,
    and copies the extracted files to the XspecT model directory.

    Args:
        url (str): The URL from which to download the models.
    """
    with TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        download_path = tmp_dir / "models.zip"
        extract_path = tmp_dir / "extracted_models"

        r = requests.get(url, allow_redirects=True, timeout=10)
        with open(download_path, "wb") as f:
            f.write(r.content)

        shutil.unpack_archive(
            download_path,
            extract_path,
            "zip",
        )

        shutil.copytree(
            extract_path,
            get_xspect_model_path(),
            dirs_exist_ok=True,
        )

        shutil.rmtree(extract_path)
