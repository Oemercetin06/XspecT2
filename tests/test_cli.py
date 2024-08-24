"""Test XspecT CLI"""

import json
import pytest
from click.testing import CliRunner
from src.xspect.main import cli


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
    result = runner.invoke(cli, ["classify", genus, assembly_file_path])

    run_path = result.output.strip().split("'")[1]
    with open(run_path, encoding="utf-8") as f:
        result_content = json.load(f)
        assert result_content["results"][0]["prediction"] == species


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
    result = runner.invoke(cli, ["classify", "-m", genus, assembly_file_path])

    run_path = result.output.strip().split("'")[1]
    with open(run_path, encoding="utf-8") as f:
        result_content = json.load(f)
        assert (
            result_content["results"][0]["subprocessing_steps"][0]["result"][
                "prediction"
            ]
            == species
        )
