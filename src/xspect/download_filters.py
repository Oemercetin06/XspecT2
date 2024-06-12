"""Download filters from public repository."""

import os
import shutil
import requests

from xspect.definitions import get_xspect_model_path, get_xspect_tmp_path


def download_test_filters(url):
    """Download filters."""

    download_path = get_xspect_tmp_path() / "models.zip"
    extract_path = get_xspect_tmp_path() / "extracted_models"

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

    os.remove(download_path)
    shutil.rmtree(extract_path)
