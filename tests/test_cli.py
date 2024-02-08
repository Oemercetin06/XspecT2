"""Test XspecT Mini (CLI)"""

import subprocess
import sys
import pytest
from click.testing import CliRunner
from src.xspect.main import cli


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
    runner = CliRunner()
    result = runner.invoke(cli, ["classify", genus, assembly_dir_path])
    assert species in result.output


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
    runner = CliRunner()
    result = runner.invoke(cli, ["classify", "-i", "Acinetobacter", assembly_dir_path])
    assert ic in result.output


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
    runner = CliRunner()
    result = runner.invoke(cli, ["classify", "-o", "Acinetobacter", assembly_dir_path])
    assert all(x in result.output for x in oxa)
