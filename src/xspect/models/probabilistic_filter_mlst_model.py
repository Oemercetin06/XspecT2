"""Probabilistic filter MLST model for sequence data"""

__author__ = "Cetin, Oemer"

import cobs_index
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cobs_index import DocumentList
from collections import defaultdict
from xspect.file_io import get_record_iterator
from xspect.mlst_feature.mlst_helper import MlstResult


class ProbabilisticFilterMlstSchemeModel:
    """Probabilistic filter MLST scheme model for sequence data"""

    def __init__(
        self,
        k: int,
        model_display_name: str,
        base_path: Path,
        fpr: float = 0.001,
    ) -> None:
        if k < 1:
            raise ValueError("Invalid k value, must be greater than 0")
        if not isinstance(base_path, Path):
            raise ValueError("Invalid base path, must be a pathlib.Path object")

        self.k = k
        self.model_display_name = model_display_name
        self.base_path = base_path / "MLST"
        self.fpr = fpr
        self.model_type = "Strain"
        self.loci = {}
        self.scheme_path = ""
        self.cobs_path = ""
        self.avg_locus_bp_size = []
        self.indices = []

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the model"""
        return {
            "k": self.k,
            "model_display_name": self.model_display_name,
            "model_type": self.model_type,
            "fpr": self.fpr,
            "scheme_path": str(self.scheme_path),
            "cobs_path": str(self.cobs_path),
            "average_locus_base_pair_size": self.avg_locus_bp_size,
            "loci": self.loci,
        }

    def get_cobs_index_path(self, scheme: str, locus: str) -> Path:
        """Returns the path to the cobs index"""
        # To differentiate from genus and species models
        cobs_path = self.base_path / f"{scheme}"
        cobs_path.mkdir(exist_ok=True, parents=True)
        return cobs_path / f"{locus}.cobs_compact"

    def fit(self, scheme_path: Path) -> None:
        """Trains a COBS structure for every locus with all its alleles"""
        if not scheme_path.exists():
            raise ValueError(
                "Scheme not found. Please make sure to download the schemes prior!"
            )

        scheme = str(scheme_path).split("/")[-1]
        cobs_path = ""
        # COBS structure for every locus (default = 7 for Oxford or Pasteur scheme)
        for locus_path in sorted(scheme_path.iterdir()):
            locus = str(locus_path).split("/")[-1]
            # counts all fasta files that belong to a locus
            self.loci[locus] = sum(
                (1 for _ in locus_path.iterdir() if not str(_).endswith("cache"))
            )

            # determine the avg base pair size of alleles
            fasta_file = next(locus_path.glob("*.fasta"), None)
            with open(fasta_file, "r") as handle:
                record = next(SeqIO.parse(handle, "fasta"))
            self.avg_locus_bp_size.append(len(record.seq))

            # COBS only accepts strings as paths
            doclist = DocumentList(str(locus_path))
            index_params = cobs_index.CompactIndexParameters()
            index_params.term_size = self.k  # k-mer size
            index_params.clobber = True  # overwrite output and temporary files
            index_params.false_positive_rate = self.fpr

            # Creates COBS data structure for each locus
            cobs_path = self.get_cobs_index_path(scheme, locus)
            cobs_index.compact_construct_list(doclist, str(cobs_path), index_params)
            # Saves COBS-file inside the "indices" attribute
            self.indices.append(cobs_index.Search(str(cobs_path)))

        self.scheme_path = scheme_path
        self.cobs_path = cobs_path.parent

    def save(self) -> None:
        """Saves the model to disk"""
        scheme = str(self.scheme_path).split("/")[
            -1
        ]  # [-1] -> contains the scheme name
        json_path = self.base_path / scheme / f"{scheme}.json"
        json_object = json.dumps(self.to_dict(), indent=4)

        with open(json_path, "w", encoding="utf-8") as file:
            file.write(json_object)

    @staticmethod
    def load(scheme_path: Path) -> "ProbabilisticFilterMlstSchemeModel":
        """Loads the model from a JSON-file"""
        scheme_name = str(scheme_path).split("/")[-1]
        json_path = scheme_path / f"{scheme_name}.json"
        with open(json_path, "r", encoding="utf-8") as file:
            json_object = file.read()
            model_json = json.loads(json_object)
            model = ProbabilisticFilterMlstSchemeModel(
                model_json["k"],
                model_json["model_display_name"],
                json_path.parent,
                model_json["fpr"],
            )
            model.scheme_path = model_json["scheme_path"]
            model.cobs_path = model_json["cobs_path"]
            model.avg_locus_bp_size = model_json["average_locus_base_pair_size"]
            model.loci = model_json["loci"]

            for entry in sorted(json_path.parent.iterdir()):
                if not entry.exists():
                    raise FileNotFoundError(f"Index file not found at {entry}")
                if str(entry).endswith(".json"):  # only COBS-files
                    continue
                model.indices.append(cobs_index.Search(str(entry), False))
            return model

    def calculate_hits(self, path: Path, sequence: Seq, step: int = 1) -> list[dict]:
        """Calculates the hits for a sequence"""
        if not isinstance(sequence, Seq):
            raise ValueError("Invalid sequence, must be a Bio.Seq object")

        if not len(sequence) > self.k:
            raise ValueError("Invalid sequence, must be longer than k")

        if not self.indices:
            raise ValueError("The Model has not been trained yet")

        scheme_path_list = []
        for entry in sorted(path.iterdir()):
            if str(entry).endswith(".json"):
                continue
            file_name = str(entry).split("/")[-1]  # file_name = locus
            scheme_path_list.append(file_name.split(".")[0])  # without the file ending

        result_dict = {}
        highest_results = {}
        counter = 0
        # split the sequence in parts based on sequence length
        if len(sequence) >= 10000:
            for index in self.indices:
                cobs_results = []
                allele_len = self.avg_locus_bp_size[counter]
                split_sequence = self.sequence_splitter(str(sequence), allele_len)
                for split in split_sequence:
                    res = index.search(split, step=step)
                    split_result = self.get_cobs_result(res)
                    if not split_result:
                        continue
                    cobs_results.append(split_result)

                all_counts = defaultdict(int)
                for result in cobs_results:
                    for name, value in result.items():
                        all_counts[name] += value

                sorted_counts = dict(
                    sorted(all_counts.items(), key=lambda item: -item[1])
                )
                first_key = next(iter(sorted_counts))
                highest_result = sorted_counts[first_key]
                result_dict[scheme_path_list[counter]] = sorted_counts
                highest_results[scheme_path_list[counter]] = {first_key: highest_result}
                counter += 1
        else:
            for index in self.indices:
                res = index.search(
                    str(sequence), step=step
                )  # COBS can't handle Seq-Objects
                result_dict[scheme_path_list[counter]] = self.get_cobs_result(res)
                highest_results[scheme_path_list[counter]] = (
                    self.get_highest_cobs_result(res)
                )
                counter += 1
        return [{"Strain type": highest_results}, {"All results": result_dict}]

    def predict(
        self,
        cobs_path: Path,
        sequence_input: (
            SeqRecord
            | list[SeqRecord]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Path
        ),
        step: int = 1,
    ) -> MlstResult:
        """Returns scores for the sequence(s) based on the filters in the model"""
        if isinstance(sequence_input, SeqRecord):
            if sequence_input.id == "<unknown id>":
                sequence_input.id = "test"
            hits = {
                sequence_input.id: self.calculate_hits(cobs_path, sequence_input.seq)
            }
            return MlstResult(self.model_display_name, step, hits)

        if isinstance(sequence_input, Path):
            return ProbabilisticFilterMlstSchemeModel.predict(
                self, cobs_path, get_record_iterator(sequence_input), step=step
            )

        if isinstance(
            sequence_input,
            (SeqIO.FastaIO.FastaIterator, SeqIO.QualityIO.FastqPhredIterator),
        ):
            hits = {}
            # individual_seq is a SeqRecord-Object
            for individual_seq in sequence_input:
                individual_hits = self.calculate_hits(cobs_path, individual_seq.seq)
                hits[individual_seq.id] = individual_hits
            return MlstResult(self.model_display_name, step, hits)

        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a"
            " SeqIO FastaIterator, or a SeqIO FastqPhredIterator"
        )

    def get_highest_cobs_result(self, cobs_result: cobs_index.SearchResult) -> dict:
        """Returns the first entry in a COBS search result."""
        # counter = 1
        # dictio = {}
        for individual_result in cobs_result:
            # COBS already sorts the result in descending order
            # The first doc_name has the highest result which is needed to determine the allele
            return {individual_result.doc_name: individual_result.score}

    def get_cobs_result(self, cobs_result: cobs_index.SearchResult) -> dict:
        """Returns all entries in a COBS search result."""
        return {
            individual_result.doc_name: individual_result.score
            for individual_result in cobs_result
            if individual_result.score > 50
        }

    def sequence_splitter(self, input_sequence: str, allele_len: int) -> list[str]:
        """Returns an equally divided sequence in form of a list."""
        # An input sequence will have 10000 or more base pairs.
        sequence_len = len(input_sequence)

        if sequence_len < 100000:
            substring_length = allele_len // 10
        elif 100000 <= sequence_len < 1000000:
            substring_length = allele_len
        elif 1000000 <= sequence_len < 10000000:
            substring_length = allele_len * 10
        else:
            substring_length = allele_len * 100

        substring_list = []
        start = 0

        while start + substring_length <= sequence_len:
            substring_list.append(input_sequence[start : start + substring_length])
            start += substring_length - self.k + 1  # To not lose kmers when dividing

        # The remaining string is either appended to the list or added to the last entry.
        if start < len(input_sequence):
            remaining_substring = input_sequence[start:]
            # A substring needs to be at least of size k for COBS.
            if len(remaining_substring) < self.k:
                substring_list[-1] += remaining_substring
            else:
                substring_list.append(remaining_substring)
        return substring_list
