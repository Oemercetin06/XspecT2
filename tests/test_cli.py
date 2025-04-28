"""Test XspecT CLI"""

import json
import pytest
from click.testing import CliRunner
from xspect.main import cli


@pytest.mark.parametrize(
    ["assembly_file_path", "genus", "species"],
    [
        (
            "GCF_000069245.1_ASM6924v1_genomic.fna",
            "Acinetobacter",
            "470",
        ),
        (
            "GCF_000018445.1_ASM1844v1_genomic.fna",
            "Acinetobacter",
            "470",
        ),
        ("GCF_000006945.2_ASM694v2_genomic.fna", "Salmonella", "28901"),
    ],
    indirect=["assembly_file_path"],
)
def test_species_assignment(assembly_file_path, genus, species):
    """Test the species assignment"""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "classify",
            "species",
            "-g",
            genus,
            "-i",
            assembly_file_path,
            "-o",
            "out.json",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"

    with open("out.json", encoding="utf-8") as f:
        result_content = json.load(f)
        assert result_content["prediction"] == species


@pytest.mark.parametrize(
    ["assembly_file_path", "genus", "species"],
    [
        (
            "GCF_000069245.1_ASM6924v1_genomic.fna",
            "Acinetobacter",
            "470",
        ),
    ],
    indirect=["assembly_file_path"],
)
def test_metagenome_mode(assembly_file_path, genus, species):
    """Test the metagenome mode"""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "filter",
            "genus",
            "-g",
            genus,
            "-i",
            assembly_file_path,
            "-o",
            "filtered.fna",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"
    result = runner.invoke(
        cli,
        [
            "classify",
            "species",
            "-g",
            genus,
            "-i",
            "filtered.fna",
            "-o",
            "out.json",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"
    with open("out.json", encoding="utf-8") as f:
        result_content = json.load(f)
        assert result_content["prediction"] == species
