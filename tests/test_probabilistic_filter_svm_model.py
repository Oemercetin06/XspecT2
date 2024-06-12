""" Tests for the ProbabilisticFilterSVMModel class. """

# pylint: disable=redefined-outer-name

from linecache import getline
from pathlib import Path
import pytest
from src.xspect.models.probabilistic_filter_svm_model import ProbabilisticFilterSVMModel


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
        f.readline()
        scores = f.readlines()
        # Check if the diagonal of the csv scores matrix is 1
        assert scores[0].split(",")[3] == "1.0"
        assert scores[1].split(",")[2] == "1.0"
        assert scores[2].split(",")[1] == "1.0"

        # Check if scores outside the diagonal are less than 1
        assert float(scores[0].split(",")[1]) < 1
        assert float(scores[0].split(",")[2]) < 1
        assert float(scores[1].split(",")[1]) < 1
        assert float(scores[1].split(",")[3]) < 1
        assert float(scores[2].split(",")[2]) < 1
        assert float(scores[2].split(",")[3]) < 1


def test_predict(trained_filter_model, multiple_assembly_svm_path):
    """Test the predict method of the ProbabilisticFilterSVMModel class."""
    for file in Path(multiple_assembly_svm_path).iterdir():
        if file.suffix not in [".fasta", ".fa", ".fna", ".fastq", ".fq"]:
            continue
        prediction, _ = trained_filter_model.predict(file)
        file_header = getline(str(file), 1)
        label_id = file_header.replace("\n", "").replace(">", "")
        assert prediction == label_id


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


def test_save_and_load(filter_model):
    """Test the load method of the ProbabilisticFilterSVMModel class."""
    filter_model.save()
    json_path = filter_model.base_path / "test-filter-species.json"
    assert (json_path).is_file()

    loaded_model = ProbabilisticFilterSVMModel.load(json_path)
    assert filter_model.to_dict() == loaded_model.to_dict()
