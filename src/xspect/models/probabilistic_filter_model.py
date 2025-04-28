"""Probabilistic filter model for sequence data"""

import json
from math import ceil
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from slugify import slugify
import cobs_index as cobs
from xspect.definitions import fasta_endings, fastq_endings
from xspect.file_io import get_record_iterator
from xspect.models.result import ModelResult


class ProbabilisticFilterModel:
    """Probabilistic filter model for sequence data"""

    def __init__(
        self,
        k: int,
        model_display_name: str,
        author: str,
        author_email: str,
        model_type: str,
        base_path: Path,
        fpr: float = 0.01,
        num_hashes: int = 7,
        training_accessions: dict[str, list[str]] = None,
    ) -> None:
        if k < 1:
            raise ValueError("Invalid k value, must be greater than 0")
        if not model_display_name:
            raise ValueError("Invalid filter display name, must be a non-empty string")
        if not model_type:
            raise ValueError("Invalid filter type, must be a non-empty string")
        if not isinstance(base_path, Path):
            raise ValueError("Invalid base path, must be a pathlib.Path object")

        self.k = k
        self.model_display_name = model_display_name
        self.author = author
        self.author_email = author_email
        self.model_type = model_type
        self.base_path = base_path
        self.display_names = {}
        self.fpr = fpr
        self.num_hashes = num_hashes
        self.index = None
        self.training_accessions = training_accessions

    def get_cobs_index_path(self) -> Path:
        """Returns the path to the cobs index"""
        return str(self.base_path / self.slug() / "index.cobs_classic")

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the model"""
        return {
            "k": self.k,
            "model_display_name": self.model_display_name,
            "author": self.author,
            "author_email": self.author_email,
            "model_type": self.model_type,
            "model_class": self.__class__.__name__,
            "display_names": self.display_names,
            "fpr": self.fpr,
            "num_hashes": self.num_hashes,
            "training_accessions": self.training_accessions,
        }

    def slug(self) -> str:
        """Returns a slug representation of the model"""
        return slugify(self.model_display_name + "-" + str(self.model_type))

    def fit(
        self,
        dir_path: Path,
        display_names: dict = None,
        training_accessions: dict[str, list[str]] = None,
    ) -> None:
        """Adds filters to the model"""

        if display_names is None:
            display_names = {}

        if not isinstance(dir_path, Path):
            raise ValueError("Invalid directory path, must be a pathlib.Path object")

        if not dir_path.exists():
            raise ValueError("Directory path does not exist")

        if not dir_path.is_dir():
            raise ValueError("Directory path must be a directory")

        self.training_accessions = training_accessions

        doclist = cobs.DocumentList()
        for file in dir_path.iterdir():
            if file.is_file() and file.suffix[1:] in fasta_endings + fastq_endings:
                # cobs only uses the file name to the first "." as the document name
                if file.stem in display_names:
                    self.display_names[file.stem.split(".")[0]] = display_names[
                        file.stem
                    ]
                else:
                    self.display_names[file.stem.split(".")[0]] = file.stem
                doclist.add(str(file))

        if len(doclist) == 0:
            raise ValueError(
                "No valid files found in directory. Must be fasta or fastq"
            )

        index_params = cobs.ClassicIndexParameters()
        index_params.term_size = self.k
        index_params.num_hashes = self.num_hashes
        index_params.false_positive_rate = self.fpr
        index_params.clobber = True

        cobs.classic_construct_list(doclist, self.get_cobs_index_path(), index_params)

        self.index = cobs.Search(self.get_cobs_index_path(), True)

    def calculate_hits(
        self, sequence: Seq, filter_ids: list[str] = None, step: int = 1
    ) -> dict:
        """Calculates the hits for a sequence"""

        if not isinstance(sequence, (Seq)):
            raise ValueError(
                "Invalid sequence, must be a Bio.Seq or a Bio.SeqRecord object"
            )

        if not len(sequence) > self.k:
            raise ValueError("Invalid sequence, must be longer than k")

        r = self.index.search(str(sequence), step=step)
        result_dict = self._convert_cobs_result_to_dict(r)
        if filter_ids:
            return {doc: result_dict[doc] for doc in filter_ids}
        return result_dict

    def predict(
        self,
        sequence_input: (
            SeqRecord
            | list[SeqRecord]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Path
        ),
        filter_ids: list[str] = None,
        step: int = 1,
    ) -> ModelResult:
        """Returns scores for the sequence(s) based on the filters in the model"""
        if isinstance(sequence_input, (SeqRecord)):
            return ProbabilisticFilterModel.predict(
                self, [sequence_input], filter_ids, step=step
            )

        if self._is_sequence_list(sequence_input) | self._is_sequence_iterator(
            sequence_input
        ):
            hits = {}
            num_kmers = {}
            for individual_sequence in sequence_input:
                individual_hits = self.calculate_hits(
                    individual_sequence.seq, filter_ids, step=step
                )
                num_kmers[individual_sequence.id] = self._count_kmers(
                    individual_sequence, step=step
                )
                hits[individual_sequence.id] = individual_hits
            return ModelResult(self.slug(), hits, num_kmers, sparse_sampling_step=step)

        if isinstance(sequence_input, Path):
            return ProbabilisticFilterModel.predict(
                self, get_record_iterator(sequence_input), step=step
            )

        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a"
            " SeqIO FastaIterator, a SeqIO FastqPhredIterator, or a Path object to a"
            " fasta/fastq file"
        )

    def save(self) -> None:
        """Saves the model to disk"""
        json_path = self.base_path / f"{self.slug()}.json"
        filter_path = self.base_path / self.slug()
        filter_path.mkdir(exist_ok=True, parents=True)

        json_object = json.dumps(self.to_dict(), indent=4)

        with open(json_path, "w", encoding="utf-8") as file:
            file.write(json_object)

    @staticmethod
    def load(path: Path) -> "ProbabilisticFilterModel":
        """Loads the model from a file"""
        with open(path, "r", encoding="utf-8") as file:
            json_object = file.read()
            model_json = json.loads(json_object)
            model = ProbabilisticFilterModel(
                model_json["k"],
                model_json["model_display_name"],
                model_json["author"],
                model_json["author_email"],
                model_json["model_type"],
                path.parent,
                model_json["fpr"],
                model_json["num_hashes"],
                model_json["training_accessions"],
            )
            model.display_names = model_json["display_names"]

            p = model.get_cobs_index_path()
            if not Path(p).exists():
                raise FileNotFoundError(f"Index file not found at {p}")
            model.index = cobs.Search(p, True)

            return model

    def _convert_cobs_result_to_dict(self, cobs_result: cobs.SearchResult) -> dict:
        return {
            individual_result.doc_name: individual_result.score
            for individual_result in cobs_result
        }

    def _count_kmers(
        self,
        sequence_input: (
            Seq
            | SeqRecord
            | list[Seq]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
        ),
        step: int = 1,
    ) -> int:
        """Counts the number of kmers in the sequence(s)"""
        if isinstance(sequence_input, Seq):
            return self._count_kmers([sequence_input], step=step)

        if isinstance(sequence_input, SeqRecord):
            return self._count_kmers(sequence_input.seq, step=step)

        is_sequence_list = isinstance(sequence_input, list) and all(
            isinstance(seq, Seq) for seq in sequence_input
        )
        is_iterator = isinstance(
            sequence_input,
            (SeqIO.FastaIO.FastaIterator, SeqIO.QualityIO.FastqPhredIterator),
        )

        if is_sequence_list | is_iterator:
            kmer_sum = 0
            for individual_sequence in sequence_input:
                # we need to look specifically at .seq for SeqIO iterators
                seq = individual_sequence.seq if is_iterator else individual_sequence
                num_kmers = ceil((len(seq) - self.k + 1) / step)
                kmer_sum += num_kmers
            return kmer_sum

        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a"
            " SeqIO FastaIterator, or a SeqIO FastqPhredIterator"
        )

    def _is_sequence_list(self, sequence_input):
        return isinstance(sequence_input, list) and all(
            isinstance(seq, (SeqRecord)) for seq in sequence_input
        )

    def _is_sequence_iterator(self, sequence_input):
        return isinstance(
            sequence_input,
            (SeqIO.FastaIO.FastaIterator, SeqIO.QualityIO.FastqPhredIterator),
        )
