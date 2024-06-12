"""
File IO module tests.
"""

import os
from pathlib import Path
from src.xspect.file_io import (
    check_folder_structure,
    concatenate_meta,
)


def test_check_folder_structure(tmpdir, monkeypatch):
    """Test if the folder structure is created correctly."""
    monkeypatch.chdir(tmpdir)

    check_folder_structure()

    path = Path(os.getcwd())
    assert os.path.isdir(path / "filter")
    assert os.path.isdir(path / "genus_metadata")
    assert os.path.isdir(path / "filter" / "array_sizes")
    assert os.path.isdir(path / "filter" / "Metagenomes")
    assert os.path.isdir(path / "filter" / "species_names")
    assert os.path.isdir(path / "filter" / "translation_dicts")


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
