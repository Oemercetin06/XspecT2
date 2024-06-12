"""Probabilistic filter SVM model for sequence data"""

# pylint: disable=no-name-in-module, too-many-instance-attributes, arguments-renamed

import csv
import json
from linecache import getline
from pathlib import Path
from sklearn.svm import SVC
from Bio.Seq import Seq
from Bio import SeqIO
from xspect.models.probabilistic_filter_model import ProbabilisticFilterModel


class ProbabilisticFilterSVMModel(ProbabilisticFilterModel):
    """Probabilistic filter SVM model for sequence data"""

    def __init__(
        self,
        k: int,
        model_display_name: str,
        author: str,
        author_email: str,
        model_type: str,
        base_path: Path,
        kernel: str,
        c: float,
        fpr: float = 0.01,
        num_hashes: int = 7,
    ) -> None:
        super().__init__(
            k=k,
            model_display_name=model_display_name,
            author=author,
            author_email=author_email,
            model_type=model_type,
            base_path=base_path,
            fpr=fpr,
            num_hashes=num_hashes,
        )
        self.kernel = kernel
        self.c = c

    def to_dict(self) -> dict:
        return super().to_dict() | {
            "kernel": self.kernel,
            "C": self.c,
        }

    def set_svm_params(self, kernel: str, c: float) -> None:
        """Set the parameters for the SVM"""
        self.kernel = kernel
        self.c = c
        self.save()

    def fit(self, dir_path: Path, svm_path: Path, display_names: dict = None) -> None:
        """Fit the SVM to the sequences and labels"""

        super().fit(dir_path, display_names=display_names)

        score_list = []
        for file in svm_path.iterdir():
            if file.suffix not in [".fasta", ".fa", ".fna", ".fastq", ".fq"]:
                continue
            print(f"Calculating {file.name} scores for SVM training...")
            scores, _ = super().predict(file)
            accession = "".join(file.name.split("_")[:2])
            file_header = getline(str(file), 1)
            label_id = file_header.replace("\n", "").replace(">", "")

            # format scores for csv
            scores = dict(sorted(scores.items()))
            scores = ",".join([str(score) for score in scores.values()])
            scores = f"{accession},{scores},{label_id}"
            score_list.append(scores)

        # csv header
        keys = list(self.display_names.keys())
        keys.sort()
        score_list.insert(0, f"file,{','.join(keys)},label_id")

        with open(
            self.base_path / self.slug() / "scores.csv", "w", encoding="utf-8"
        ) as file:
            file.write("\n".join(score_list))

    def predict(
        self,
        sequence_input: (
            Seq
            | list[Seq]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Path
        ),
        filter_ids: list[str] = None,
    ) -> dict:
        """Predict the labels of the sequences"""
        # get scores and format them for the SVM
        scores, _ = super().predict(sequence_input, filter_ids)
        svm_scores = dict(sorted(scores.items()))
        svm_scores = [list(svm_scores.values())]

        svm = self._get_svm(filter_ids)
        return svm.predict(svm_scores), scores

    def _get_svm(self, id_keys) -> SVC:
        """Get the SVM for the given id keys"""
        svm = SVC(kernel=self.kernel, C=self.c)
        # parse csv
        with open(
            self.base_path / self.slug() / "scores.csv", "r", encoding="utf-8"
        ) as file:
            file.readline()
            x_train = []
            y_train = []
            for row in csv.reader(file):
                if id_keys is None or row[-1] in id_keys:
                    x_train.append(row[1:-1])
                    y_train.append(row[-1])

        # train svm
        svm.fit(x_train, y_train)
        return svm

    @staticmethod
    def load(path: Path) -> "ProbabilisticFilterSVMModel":
        """Load the model from disk"""
        with open(path, "r", encoding="utf-8") as file:
            json_object = file.read()
            model_json = json.loads(json_object)
            model = ProbabilisticFilterSVMModel(
                model_json["k"],
                model_json["model_display_name"],
                model_json["author"],
                model_json["author_email"],
                model_json["model_type"],
                path.parent,
                model_json["kernel"],
                model_json["C"],
                fpr=model_json["fpr"],
                num_hashes=model_json["num_hashes"],
            )
            model.display_names = model_json["display_names"]
            return model
