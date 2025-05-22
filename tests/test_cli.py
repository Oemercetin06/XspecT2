"""Test XspecT CLI"""

import json
from pathlib import Path
import pytest
from click.testing import CliRunner
from xspect.main import cli


def test_list_models():
    """Test the list models command"""
    runner = CliRunner()
    result = runner.invoke(cli, ["models", "list"])
    assert result.exit_code == 0, f"Error: {result.output}"
    assert "Genus" in result.output
    assert "Species" in result.output


@pytest.mark.parametrize(
    "assembly_file_path",
    [
        "GCF_000069245.1_ASM6924v1_genomic.fna",
    ],
    indirect=["assembly_file_path"],
)
def test_classify_genus(assembly_file_path, tmpdir):
    """Test the classify genus command"""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "classify",
            "genus",
            "-g",
            "Acinetobacter",
            "-i",
            assembly_file_path,
            "-o",
            str(tmpdir) + "/classify_genus.json",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"
    with open(str(tmpdir) + "/classify_genus.json", encoding="utf-8") as f:
        result_content = json.load(f)
        assert result_content["scores"]["total"]["Acinetobacter"] == 0.85


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
def test_classify_species(assembly_file_path, genus, species, tmpdir):
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
            str(tmpdir) + "/classify_species.json",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"

    with open(str(tmpdir) + "/classify_species.json", encoding="utf-8") as f:
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
def test_filter_genus_and_classify_species(assembly_file_path, genus, species, tmpdir):
    """Test filtering by a genus and then classifying species ("metagenome mode")"""
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
            str(tmpdir) + "/genus_filtered.fna",
            "-t",
            "0.7",
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
            str(tmpdir) + "/genus_filtered.fna",
            "-o",
            str(tmpdir) + "/out.json",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"
    with open(str(tmpdir) + "/out.json", encoding="utf-8") as f:
        result_content = json.load(f)
        assert result_content["prediction"] == species


@pytest.mark.parametrize(
    "assembly_file_path",
    [
        "GCF_000006945.2_ASM694v2_genomic.fna",
    ],
    indirect=["assembly_file_path"],
)
def test_filter_species(assembly_file_path, tmpdir):
    """Test filtering by species"""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "filter",
            "species",
            "-g",
            "Salmonella",
            "-s",
            "enterica",
            "-i",
            assembly_file_path,
            "-o",
            str(tmpdir) + "/species_filtered.fna",
            "-t",
            "0.7",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"
    assert Path(str(tmpdir) + "/species_filtered.fna").exists()


def test_filter_species_max_scoring(mixed_species_assembly_file_path, tmpdir):
    """Test filtering by species"""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "filter",
            "species",
            "-g",
            "Acinetobacter",
            "-s",
            "calcoaceticus",
            "-i",
            mixed_species_assembly_file_path,
            "-o",
            str(tmpdir) + "/mixed_species_filtered.fna",
            "-t",
            "-1",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"
    assert Path(str(tmpdir) + "/mixed_species_filtered.fna").exists()

    with open(str(tmpdir) + "/mixed_species_filtered.fna", encoding="utf-8") as f:
        filtered_content = f.read()
        assert "Acinetobacter calcoaceticus" in filtered_content
        assert "Acinetobacter baumannii" not in filtered_content
