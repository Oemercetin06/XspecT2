"""
File IO module tests.
"""

from pathlib import Path
from xspect.file_io import concatenate_metagenome


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
        "test1.fasta",
        "test2.fasta",
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
    meta_file = Path(tmpdir) / "Test.fasta"

    # Call the function to be tested
    concatenate_metagenome(concatenate_dir, meta_file)

    with open(meta_file, "r", encoding="utf-8") as f:
        content = f.read()
        assert content.startswith(">")
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
