"""Download filters from public repository."""

import os
import shutil
import requests


def download_test_filters(url):
    """Download filters."""

    if not os.path.exists("filter"):
        os.makedirs("filter")

    if not os.path.exists("Training_data"):
        os.makedirs("Training_data")

    r = requests.get(url, allow_redirects=True, timeout=10)
    with open("filter/filters.zip", "wb") as f:
        f.write(r.content)

    shutil.unpack_archive(
        "filter/filters.zip",
        "filter/temp",
        "zip",
    )

    shutil.copytree(
        "filter/temp/filters/Training_data",
        "Training_data",
        dirs_exist_ok=True,
    )

    shutil.rmtree("filter/temp/filters/Training_data")

    shutil.copytree(
        "filter/temp/filters",
        "filter",
        dirs_exist_ok=True,
    )

    shutil.rmtree("filter/temp")

    os.remove("filter/filters.zip")
