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
    ["assembly_file_path", "genus", "species_display_name"],
    [
        (
            "GCF_000069245.1_ASM6924v1_genomic.fna",
            "Acinetobacter",
            "470 - baumannii",
        ),
    ],
    indirect=["assembly_file_path"],
)
def test_classify_species_with_names(
    assembly_file_path, genus, species_display_name, tmpdir
):
    """Test the species assignment with display names included"""
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
            str(tmpdir) + "/classify_species_with_names.json",
            "-n",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"

    with open(str(tmpdir) + "/classify_species_with_names.json", encoding="utf-8") as f:
        result_content = json.load(f)
        assert (
            result_content["prediction"] == species_display_name.split("-")[0].strip()
        )
        for subseq_scores in result_content["scores"].values():
            assert species_display_name in subseq_scores
        for subseq_hits in result_content["hits"].values():
            assert species_display_name in subseq_hits


@pytest.mark.parametrize(
    "assembly_file_path",
    [
        "GCF_000069245.1_ASM6924v1_genomic.fna",
    ],
    indirect=["assembly_file_path"],
)
def test_classify_species_with_exclusion(assembly_file_path, tmpdir):
    """Test the species assignment with species exclusion"""
    runner = CliRunner()
    excluded_species_id = "48296"
    result = runner.invoke(
        cli,
        [
            "classify",
            "species",
            "-g",
            "Acinetobacter",
            "-i",
            assembly_file_path,
            "-o",
            str(tmpdir) + "/classify_species_with_exclusion.json",
            "--exclude-species",
            excluded_species_id,
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"

    with open(
        str(tmpdir) + "/classify_species_with_exclusion.json", encoding="utf-8"
    ) as f:
        result_content = json.load(f)
        assert result_content["prediction"] != excluded_species_id
        for subseq_scores in result_content["scores"].values():
            assert excluded_species_id not in subseq_scores
        for subseq_hits in result_content["hits"].values():
            assert excluded_species_id not in subseq_hits


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
def test_all_pipeline(assembly_file_path, genus, species, tmpdir):
    """Test the complete 'xspect all' pipeline"""
    runner = CliRunner()
    output_dir = str(tmpdir) + "/all_results"

    result = runner.invoke(
        cli,
        [
            "all",
            "-g",
            genus,
            "-i",
            assembly_file_path,
            "-o",
            output_dir,
            "-t",
            "0.7",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"

    output_path = Path(output_dir)
    assert output_path.exists()

    filtered_dir = output_path / "filtered_sequences"
    assert filtered_dir.exists()

    filtered_files = list(filtered_dir.glob("genus_filtered_*.fasta"))
    assert len(filtered_files) > 0, "No genus filtered files found"

    genus_classification_files = list(output_path.glob("genus_classification_*.json"))
    assert len(genus_classification_files) > 0, "No genus classification files found"

    species_classification_files = list(
        output_path.glob("species_classification_*.json")
    )
    assert (
        len(species_classification_files) > 0
    ), "No species classification files found"

    with open(species_classification_files[0], encoding="utf-8") as f:
        result_content = json.load(f)
        assert result_content["prediction"] == species

    if species == "470":
        mlst_classification_files = list(output_path.glob("mlst_classification_*.json"))
        assert (
            len(mlst_classification_files) > 0
        ), "No MLST classification files found for abaumannii"


@pytest.mark.parametrize(
    ["assembly_file_path", "genus", "species"],
    [
        (
            "GCF_000006945.2_ASM694v2_genomic.fna",
            "Salmonella",
            "28901",
        ),
    ],
    indirect=["assembly_file_path"],
)
def test_all_pipeline_no_mlst(assembly_file_path, genus, species, tmpdir):
    """Test 'xspect all' pipeline without MLST (non-abaumannii species)"""
    runner = CliRunner()
    output_dir = str(tmpdir) + "/all_results_no_mlst"

    result = runner.invoke(
        cli,
        [
            "all",
            "-g",
            genus,
            "-i",
            assembly_file_path,
            "-o",
            output_dir,
            "-t",
            "0.7",
        ],
    )
    assert result.exit_code == 0, f"Error: {result.output}"

    # Verify output directory structure
    output_path = Path(output_dir)
    assert output_path.exists()

    # Check for species classification results
    species_classification_files = list(
        output_path.glob("species_classification_*.json")
    )
    assert len(species_classification_files) > 0

    # Verify species classification content
    with open(species_classification_files[0], encoding="utf-8") as f:
        result_content = json.load(f)
        assert result_content["prediction"] == species

    # No MLST classification should exist for non-abaumannii species
    mlst_classification_files = list(output_path.glob("mlst_classification_*.json"))
    assert (
        len(mlst_classification_files) == 0
    ), "MLST classification should not exist for non-abaumannii"


@pytest.mark.parametrize(
    "assembly_file_path",
    [
        "GCF_000069245.1_ASM6924v1_genomic.fna",
    ],
    indirect=["assembly_file_path"],
)
def test_all_pipeline_auto_directory(assembly_file_path, tmpdir):
    """Test 'xspect all' with auto-generated output directory"""
    runner = CliRunner()

    # Change to tmpdir to ensure auto-generated directory is created there
    import os

    original_dir = os.getcwd()
    os.chdir(tmpdir)

    try:
        result = runner.invoke(
            cli,
            [
                "all",
                "-g",
                "Acinetobacter",
                "-i",
                assembly_file_path,
                "-t",
                "0.7",
            ],
        )
        assert result.exit_code == 0, f"Error: {result.output}"

        # Check that an auto-generated directory was created
        result_dirs = [
            d
            for d in Path(tmpdir).iterdir()
            if d.is_dir() and d.name.startswith("xspect_results_")
        ]
        assert len(result_dirs) > 0, "No auto-generated output directory found"

        # Verify the directory contains expected files
        output_path = result_dirs[0]
        assert (output_path / "filtered_sequences").exists()
        assert len(list(output_path.glob("genus_classification_*.json"))) > 0
        assert len(list(output_path.glob("species_classification_*.json"))) > 0
    finally:
        os.chdir(original_dir)
