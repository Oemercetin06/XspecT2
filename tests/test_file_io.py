"""
File IO module tests.
"""

import os
from pathlib import Path
import pytest
from src.xspect.file_io import (
    check_folder_structure,
    delete_non_fasta,
    get_accessions,
    get_file_paths,
    concatenate_meta,
)


def test_check_folder_structure(tmpdir):
    """Test if the folder structure is created correctly."""
    # Set up temporary directory
    os.chdir(tmpdir)

    # Call the function to be tested
    check_folder_structure()

    # Check if the folders have been created
    path = Path(os.getcwd())
    assert os.path.isdir(path / "filter")
    assert os.path.isdir(path / "genus_metadata")
    assert os.path.isdir(path / "filter" / "array_sizes")
    assert os.path.isdir(path / "filter" / "Metagenomes")
    assert os.path.isdir(path / "filter" / "species_names")
    assert os.path.isdir(path / "filter" / "translation_dicts")


@pytest.mark.parametrize(
    ["files", "expected_result"],
    [
        ([], []),
        (
            ["file1.fasta", "file2.fna", "file3.fa", "file4.ffn", "file5.frn"],
            ["file1.fasta", "file2.fna", "file3.fa", "file4.ffn", "file5.frn"],
        ),
        (
            ["file1.txt", "file2.jpg", "file3.png"],
            [],
        ),
        (
            [
                "file1.fasta",
                "file2.fna",
                "file3.fa",
                "file4.ffn",
                "file5.frn",
                "file6.txt",
                "file7.jpg",
                "file8.png",
            ],
            ["file1.fasta", "file2.fna", "file3.fa", "file4.ffn", "file5.frn"],
        ),
    ],
)
def test_delete_non_fasta(files, expected_result):
    """Test if the function deletes non-fasta files correctly."""
    assert delete_non_fasta(files) == expected_result


@pytest.mark.parametrize(
    "files, expected_result",
    [
        ([], []),
        (
            [
                "GCF_000439255.1_ASM43925v1_genomic.fna",
                "GCF_000439255.1_ASM43925v1_genomic.fna",
            ],
            ["GCF_000439255.1", "GCF_000439255.1"],
        ),
    ],
)
def test_get_accessions(files, expected_result):
    """Test if the function extracts accessions correctly."""
    assert get_accessions(files) == expected_result


@pytest.mark.parametrize(
    ["base_path", "file_names", "expected_paths"],
    [
        (Path("test"), [], []),
        (
            Path("test"),
            ["file1", "file2", "file3"],
            [Path("test/file1"), Path("test/file2"), Path("test/file3")],
        ),
    ],
)
def test_get_file_paths(base_path, file_names, expected_paths):
    """Test if the function returns the correct file paths."""
    assert get_file_paths(base_path, file_names) == expected_paths


def test_concatenate_meta(tmpdir):
    """Test if the function concatenates fasta files correctly."""
    # Set up temporary directory
    os.chdir(tmpdir)

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
        with open(concatenate_dir / file, "w") as f:
            f.write(f">{file}\n{file}")

    # Call the function to be tested
    concatenate_meta(tmpdir, "Salmonella")

    # Check if the meta file has been created and contains the correct content
    meta_file = Path(tmpdir) / "Salmonella.fasta"
    assert meta_file.is_file()

    with open(meta_file, "r") as f:
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
