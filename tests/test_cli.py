"""Test XspecT Mini (CLI)"""

import pytest
from click.testing import CliRunner
from src.xspect.main import cli


@pytest.mark.parametrize(
    ["assembly_file_path", "genus", "species"],
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
    indirect=["assembly_file_path"],
)
def test_species_assignment(assembly_file_path, genus, species):
    """Test the species assignment"""
    runner = CliRunner()
    result = runner.invoke(cli, ["classify", genus, assembly_file_path])
    assert species in result.output


@pytest.mark.parametrize(
    ["assembly_file_path", "genus", "species"],
    [
        (
            "GCF_000069245.1_ASM6924v1_genomic.fna",
            "Acinetobacter",
            "Acinetobacter baumannii",
        ),
    ],
    indirect=["assembly_file_path"],
)
def test_metagenome_mode(assembly_file_path, genus, species):
    """Test the metagenome mode"""
    runner = CliRunner()
    result = runner.invoke(cli, ["classify", "-m", genus, assembly_file_path])
    assert species in result.output
