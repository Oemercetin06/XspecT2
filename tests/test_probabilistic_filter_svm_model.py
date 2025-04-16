"""Tests for the ProbabilisticFilterSVMModel class."""

# pylint: disable=redefined-outer-name

from pathlib import Path
from csv import DictReader
import pytest
from xspect.models.probabilistic_filter_svm_model import ProbabilisticFilterSVMModel


@pytest.fixture
def filter_model(tmpdir):
    """Fixture for the ProbabilisticFilterSVMModel class."""
    base_path = Path(tmpdir.mkdir("xspect_data"))
    return ProbabilisticFilterSVMModel(
        k=21,
        model_display_name="Test Filter",
        author="John Doe",
        author_email="john.doe@example.com",
        model_type="Species",
        base_path=base_path,
        kernel="linear",
        c=1.0,
    )


@pytest.fixture
def trained_filter_model(
    filter_model, multiple_assembly_dir_path, multiple_assembly_svm_path
):
    """Fixture for the ProbabilisticFilterSVMModel class with trained model."""
    filter_model.fit(Path(multiple_assembly_dir_path), Path(multiple_assembly_svm_path))
    return filter_model


def test_fit(filter_model, multiple_assembly_dir_path, multiple_assembly_svm_path):
    """Test the fit method of the ProbabilisticFilterSVMModel class."""
    filter_model.fit(Path(multiple_assembly_dir_path), Path(multiple_assembly_svm_path))

    # Check if the scores.csv file has been created
    scores_file = filter_model.base_path / "test-filter-species" / "scores.csv"
    assert scores_file.is_file()

    with open(scores_file, "r", encoding="utf-8") as f:
        reader = DictReader(f)

        num_rows = 0

        for row in reader:
            num_rows += 1
            label_id = row["label_id"]
            assert row[label_id] == "1.0"

            # scores outside of the diagonal should be less than 1
            for key, value in row.items():
                if key not in ["label_id", "file", label_id]:
                    assert float(value) < 1

        assert num_rows == 3


def test_predict(trained_filter_model, multiple_assembly_svm_path):
    """Test the predict method of the ProbabilisticFilterSVMModel class."""
    for file in Path(multiple_assembly_svm_path).glob("**/*"):
        if file.suffix not in [".fasta", ".fa", ".fna", ".fastq", ".fq"]:
            continue
        res = trained_filter_model.predict(file, step=500)
        prediction = res.prediction
        assert prediction == file.parent.name


def test_set_svm_params(filter_model):
    """Test the set_svm_params method of the ProbabilisticFilterSVMModel class."""
    filter_model.set_svm_params(kernel="rbf", c=0.5)
    assert filter_model.kernel == "rbf"
    assert filter_model.c == 0.5


def test_save(filter_model):
    """Test the save method of the ProbabilisticFilterSVMModel class."""
    filter_model.save()
    json_path = filter_model.base_path / "test-filter-species.json"
    assert (json_path).is_file()
    with open(json_path, "r", encoding="utf-8") as f:
        data = f.read()
        assert "linear" in data
        assert "1.0" in data


def test_save_and_load(trained_filter_model):
    """Test the load method of the ProbabilisticFilterSVMModel class."""
    trained_filter_model.save()
    json_path = trained_filter_model.base_path / "test-filter-species.json"
    assert (json_path).is_file()

    loaded_model = ProbabilisticFilterSVMModel.load(json_path)
    assert trained_filter_model.to_dict() == loaded_model.to_dict()
