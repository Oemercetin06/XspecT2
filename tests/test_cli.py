"""Test XspecT Mini (CLI)"""

import subprocess
import sys
import pytest


@pytest.mark.parametrize(
    ["assembly_dir_path", "genus", "species"],
    [
        (
            "GCF_000069245.1_ASM6924v1_genomic.fna",
            "Acinetobacter",
            "Acinetobacter baumannii",
        ),
        (
            "GCF_000018445.1_ASM1844v1_genomic.fna",
            "Acinetobacter",
            "Acinetobacter baumannii",
        ),
        ("GCF_000006945.2_ASM694v2_genomic.fna", "Salmonella", "Salmonella enterica"),
    ],
    indirect=["assembly_dir_path"],
)
def test_species_assignment(assembly_dir_path, genus, species):
    """Test the species assignment"""
    result = subprocess.run(
        sys.executable
        + " XspecT_mini.py "
        + genus
        + " xspect fna "
        + assembly_dir_path,
        shell=True,
        capture_output=True,
        check=False,
    )
    assert species in result.stdout.decode("utf-8")


@pytest.mark.parametrize(
    ["assembly_dir_path", "ic"],
    [
        ("GCF_000069245.1_ASM6924v1_genomic.fna", "IC1"),
        ("GCF_000018445.1_ASM1844v1_genomic.fna", "IC2"),
    ],
    indirect=["assembly_dir_path"],
)
def test_ic_assignment(assembly_dir_path, ic):
    """Test the international clonal (IC) type assignment"""
    path = "tests/test_assemblies"
    result = subprocess.run(
        sys.executable
        + " XspecT_mini.py Acinetobacter classt fna "
        + assembly_dir_path,
        shell=True,
        capture_output=True,
        check=False,
    )
    assert ic in result.stdout.decode("utf-8")


@pytest.mark.parametrize(
    ["assembly_dir_path", "oxa"],
    [
        ("GCF_000069245.1_ASM6924v1_genomic.fna", ["OXA-51_family", "OXA-69"]),
        (
            "GCF_000018445.1_ASM1844v1_genomic.fna",
            ["OXA-51_family", "OXA-58_family", "OXA-66"],
        ),
    ],
    indirect=["assembly_dir_path"],
)
def test_oxa_assignment(assembly_dir_path, oxa):
    """Test the OXA type assignment"""
    path = "tests/test_assemblies"
    result = subprocess.run(
        sys.executable + " XspecT_mini.py Acinetobacter oxa fna " + assembly_dir_path,
        shell=True,
        capture_output=True,
        check=False,
    )
    assert all(x in result.stdout.decode("utf-8") for x in oxa)
