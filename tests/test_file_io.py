"""
File IO module tests.
"""

import os
from pathlib import Path
from xspect.file_io import (
    concatenate_meta,
)


def test_concatenate_meta(tmpdir, monkeypatch):
    """Test if the function concatenates fasta files correctly."""
    # Set up temporary directory
    monkeypatch.chdir(tmpdir)

    # Create a temporary directory for the concatenated fasta files
    concatenate_dir = Path(tmpdir) / "concatenate"
    concatenate_dir.mkdir()

    # Create some temporary fasta files
    fasta_files = [
        "file1.fasta",
        "file2.fna",
        "file3.fa",
        "file4.ffn",
        "file5.frn",
        "file6.txt",
        "file7.jpg",
        "file8.png",
    ]
    for file in fasta_files:
        with open(concatenate_dir / file, "w", encoding="utf-8") as f:
            f.write(f">{file}\n{file}")

    # Call the function to be tested
    concatenate_meta(tmpdir, "Salmonella")

    # Check if the meta file has been created and contains the correct content
    meta_file = Path(tmpdir) / "Salmonella.fasta"
    assert meta_file.is_file()

    with open(meta_file, "r", encoding="utf-8") as f:
        content = f.read()
        assert content.startswith(">Salmonella metagenome")
        for file in fasta_files:
            if (
                file.endswith(".fasta")
                or file.endswith(".fna")
                or file.endswith(".fa")
                or file.endswith(".ffn")
                or file.endswith(".frn")
            ):
                assert file in content
            else:
                assert file not in content
