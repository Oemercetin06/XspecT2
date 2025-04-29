"""Pytest fixtures for the tests in the tests/ directory"""

# pylint: disable=line-too-long

import shutil
import os
from pathlib import Path
import requests
import pytest


def pytest_sessionstart():
    """Download assemblies from NCBI"""
    assemblies = {
        "GCF_000006945.2_ASM694v2_genomic.fna": "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000006945.2/download?include_annotation_type=GENOME_FASTA",
        "GCF_000018445.1_ASM1844v1_genomic.fna": "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000018445.1/download?include_annotation_type=GENOME_FASTA",
        "GCF_000069245.1_ASM6924v1_genomic.fna": "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000069245.1/download?include_annotation_type=GENOME_FASTA",
        "GCA_900444805.1_58932_B01_genomic.fna": "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_900444805.1/download?include_annotation_type=GENOME_FASTA",
    }
    if not os.path.exists("tests/test_assemblies"):
        os.makedirs("tests/test_assemblies")

    for assembly, url in assemblies.items():
        if not os.path.exists("tests/test_assemblies/" + assembly):
            print("Downloading " + assembly)

            r = requests.get(url, allow_redirects=True, timeout=10)
            with open("tests/test_assemblies/" + assembly + ".zip", "wb") as f:
                f.write(r.content)

            # Unzip and move
            shutil.unpack_archive(
                "tests/test_assemblies/" + assembly + ".zip",
                "tests/test_assemblies/temp",
                "zip",
            )
            refseq_id = "_".join(assembly.split("_", 2)[:2])
            shutil.move(
                "tests/test_assemblies/temp/ncbi_dataset/data/"
                + refseq_id
                + "/"
                + assembly,
                "tests/test_assemblies/" + assembly,
            )

            # Clean up
            shutil.rmtree("tests/test_assemblies/temp")
            os.remove("tests/test_assemblies/" + assembly + ".zip")


@pytest.fixture
def assembly_dir_path(tmp_path, request):
    """Create a temporary directory with requested test assembly and return the path as string"""
    shutil.copy("tests/test_assemblies/" + request.param, tmp_path)
    return tmp_path.as_posix()


@pytest.fixture
def assembly_file_path(tmp_path, request):
    """Create a temporary directory with requested test assembly and returns the path to the file"""
    file_path = shutil.copy("tests/test_assemblies/" + request.param, tmp_path)
    return file_path


@pytest.fixture
def multiple_assembly_dir_path(tmp_path):
    """Create a temporary directory with multiple test assemblies and return the path as string"""
    assemblies = [
        "GCF_000006945.2_ASM694v2_genomic.fna",
        "GCF_000018445.1_ASM1844v1_genomic.fna",
        "GCF_000069245.1_ASM6924v1_genomic.fna",
    ]
    for assembly in assemblies:
        shutil.copy("tests/test_assemblies/" + assembly, tmp_path)
    return tmp_path.as_posix()


@pytest.fixture
def multiple_assembly_svm_path(tmp_path):
    """Create a temporary directory with multiple test assemblies and return the path as string"""

    svm_path = Path(tmp_path) / "svm"
    svm_path.mkdir()
    assemblies = [
        "GCF_000006945.2_ASM694v2_genomic.fna",
        "GCF_000018445.1_ASM1844v1_genomic.fna",
        "GCF_000069245.1_ASM6924v1_genomic.fna",
    ]
    for assembly in assemblies:
        output_path = svm_path / assembly.split(".", maxsplit=1)[0] / assembly
        output_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy("tests/test_assemblies/" + assembly, output_path)
    return svm_path.as_posix()


@pytest.fixture
def concatenated_assembly_file_path(tmp_path):
    """Create a temporary directory with multiple test assemblies and return the path as string"""
    # two acinetobacter assemblies
    assemblies = [
        "GCF_000018445.1_ASM1844v1_genomic.fna",
        "GCF_000069245.1_ASM6924v1_genomic.fna",
    ]
    with open(tmp_path / "concatenated_assembly.fna", "w", encoding="utf-8") as outfile:
        for assembly in assemblies:
            with open(
                "tests/test_assemblies/" + assembly, "r", encoding="utf-8"
            ) as infile:
                shutil.copyfileobj(infile, outfile)
    return (tmp_path / "concatenated_assembly.fna").as_posix()


@pytest.fixture
def mixed_species_assembly_file_path(tmp_path):
    """Create a temporary directory a fasta file which contains two mixed species assemblies"""
    # two acinetobacter assemblies
    assemblies = [
        "GCF_000018445.1_ASM1844v1_genomic.fna",
        "GCA_900444805.1_58932_B01_genomic.fna",
    ]
    with open(
        tmp_path / "mixed_species_assembly.fna", "w", encoding="utf-8"
    ) as outfile:
        for assembly in assemblies:
            with open(
                "tests/test_assemblies/" + assembly, "r", encoding="utf-8"
            ) as infile:
                shutil.copyfileobj(infile, outfile)
    return (tmp_path / "mixed_species_assembly.fna").as_posix()
