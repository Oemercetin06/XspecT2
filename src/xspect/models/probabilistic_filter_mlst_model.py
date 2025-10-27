"""Probabilistic filter MLST model for sequence data"""

# pylint: disable=arguments-differ

__author__ = "Cetin, Oemer"

import json
from pathlib import Path
from collections import defaultdict
import cobs_index
from cobs_index import DocumentList
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from slugify import slugify
from xspect.file_io import get_record_iterator
from xspect.models.mlst_result import MlstResult
from xspect.handlers.pubmlst import PubMLSTHandler
from xspect.models.probabilistic_filter_model import ProbabilisticFilterModel


class ProbabilisticFilterMlstSchemeModel(ProbabilisticFilterModel):
    """Probabilistic filter MLST scheme model for sequence data"""

    def __init__(
        self,
        k: int,
        model_display_name: str,
        base_path: Path,
        scheme_url: str,
        organism: str,
        fpr: float = 0.001,
        num_hashes: int = 1,
        author: str | None = None,
        author_email: str | None = None,
        model_type: str = "MLST",
    ) -> None:
        """Initialize a ProbabilisticFilterMlstSchemeModel object."""
        super().__init__(
            k,
            model_display_name,
            author,
            author_email,
            model_type,
            base_path,
            fpr,
            num_hashes,
            None,
        )
        self.organism = organism
        self.scheme_url = scheme_url

        self.loci = {}
        self.avg_locus_bp_size = []
        self.indices = []

    def to_dict(self) -> dict:
        """
        Convert the model to a dictionary representation

        Returns:
            dict: A dictionary containing the model's parameters and state.
        """
        return super().to_dict() | {
            "organism": self.organism,
            "scheme_url": self.scheme_url,
            "loci": self.loci,
            "average_locus_base_pair_size": self.avg_locus_bp_size,
        }

    def slug(self) -> str:
        """
        Returns a slug representation of the model

        This method generates a slug based on the model's display name and type,
        which can be used for file naming or identification purposes.

        Returns:
            str: A slug representation of the model.
        """
        return slugify(
            self.organism + "-" + self.model_display_name + "-" + self.model_type
        )

    def get_cobs_index_path(self, locus: str) -> Path:
        """
        Get the path to the cobs indices.

        This function creates a directory based on the scheme name, if it does not exist.
        A COBS-Index file is created for every locus in a scheme.

        Args:
            locus (str): The name of the locus.

        Returns:
            Path: The path to the COBS indices.
        """
        cobs_path = self.base_path / self.slug()
        return cobs_path / f"{locus}.cobs_compact"

    def fit(self, scheme_path: Path) -> None:
        """
        Trains a COBS index for every locus with all its alleles.

        This function creates COBS-indices.
        Many attributes of an object are set in this function.

        Args:
            scheme_path (Path): The path to the scheme directory with all loci.

        Raises:
            ValueError: If the scheme alleles have not been downloaded prior.
        """
        if not scheme_path.exists():
            raise ValueError(
                "Scheme not found. Please make sure to download the schemes prior!"
            )

        for locus_path in sorted(scheme_path.iterdir()):
            locus = locus_path.name
            # counts all fasta files that belong to a locus
            self.loci[locus] = sum(
                (1 for p in locus_path.iterdir() if p.suffix == ".fasta")
            )

            # determine the avg base pair size of alleles
            fasta_file_path = next(locus_path.glob("*.fasta"), None)
            with open(fasta_file_path, "r", encoding="utf-8") as fasta_file:
                record = next(SeqIO.parse(fasta_file, "fasta"))
            self.avg_locus_bp_size.append(len(record.seq))

            doclist = DocumentList(str(locus_path))
            index_params = cobs_index.CompactIndexParameters()
            index_params.term_size = self.k  # k-mer size
            index_params.clobber = True  # overwrite output and temporary files
            index_params.false_positive_rate = self.fpr
            index_params.num_hashes = self.num_hashes

            cobs_path = self.get_cobs_index_path(locus)
            cobs_path.mkdir(exist_ok=True, parents=True)
            cobs_index.compact_construct_list(doclist, str(cobs_path), index_params)
            self.indices.append(cobs_index.Search(str(cobs_path)))

    def save(self) -> None:
        """Saves the model to disk"""
        json_path = self.base_path / f"{self.slug()}.json"
        json_object = json.dumps(self.to_dict(), indent=4)

        with open(json_path, "w", encoding="utf-8") as file:
            file.write(json_object)

    @staticmethod
    def load(path: Path) -> "ProbabilisticFilterMlstSchemeModel":
        """
        Loads the model from a JSON-file.

        Args:
            path (Path): The path of the MLST model.

        Returns:
            ProbabilisticFilterMlstSchemeModel: A model from the disk in JSON format.
        """
        if not path.exists():
            raise FileNotFoundError(f"Model JSON not found at {path}")

        with open(path, "r", encoding="utf-8") as file:
            model_json = json.loads(file.read())

            model = ProbabilisticFilterMlstSchemeModel(
                model_json["k"],
                model_json["model_display_name"],
                path.parent,
                model_json["scheme_url"],
                model_json["organism"],
                model_json["fpr"],
                model_json["num_hashes"],
                model_json.get("author"),
                model_json.get("author_email"),
                model_json.get("model_type"),
            )
            model.avg_locus_bp_size = model_json.get("average_locus_base_pair_size", [])
            model.loci = model_json.get("loci", {})

            for locus in model.loci.keys():
                index_path = model.get_cobs_index_path(locus)
                if not index_path.exists():
                    raise FileNotFoundError(f"Index file not found at {index_path}")
                model.indices.append(cobs_index.Search(str(index_path), False))

            return model

    def calculate_hits(
        self,
        sequence: Seq,
        step: int = 1,
        limit: bool = False,
        limit_number: int = 5,
    ) -> list[dict]:
        """
        Calculates the hits for a sequence.

        This function has two ways of identifying strain types.
        Sequences with a length of up to 10000 base pairs are handled without preprocessing.
        Sequences with a length >= 10000 base pairs are divided into substrings.
        The results of each substring are added up to find the strain type.

        Args:
            sequence (Seq): The input sequence for classification.
            step (int, optional): The amount of kmers that are passed; defaults to one.
            limit (bool): Applying a filter that limits the best result.
            limit_number (int): The amount of results when the filter is set to true.

        Returns:
            list[dict]: The results of the prediction.

        Raises:
            ValueError: If the model has not been trained.
            ValueError: If the sequence is shorter than k.
            ValueError: If the sequence is not a Seq-object.
        """
        if not isinstance(sequence, Seq):
            raise ValueError("Invalid sequence, must be a Bio.Seq object")

        if not len(sequence) > self.k:
            raise ValueError("Invalid sequence, must be longer than k")

        if not self.indices:
            raise ValueError("The model has not been trained yet")

        scheme_path_list = list(self.loci.keys())

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
                    split_result = self.get_cobs_result(res, True)
                    if not split_result:
                        continue
                    cobs_results.append(split_result)

                # add all split results of an Allele id into one
                all_counts = defaultdict(int)
                for result in cobs_results:
                    for name, value in result.items():
                        all_counts[name] += value

                sorted_counts = dict(
                    sorted(all_counts.items(), key=lambda item: -item[1])
                )

                if limit:
                    sorted_counts = dict(list(sorted_counts.items())[:limit_number])

                if not sorted_counts:
                    result_dict = "A Strain type could not be detected because of no kmer matches!"
                    highest_results[scheme_path_list[counter]] = {"N/A": 0}
                else:
                    first_key = next(iter(sorted_counts))
                    highest_result = sorted_counts[first_key]
                    result_dict[scheme_path_list[counter]] = sorted_counts
                    highest_results[scheme_path_list[counter]] = {
                        first_key: highest_result
                    }
                counter += 1
        else:  # No split procedure is needed, when the sequence is short
            for index in self.indices:
                res = index.search(  # COBS can't handle Seq-Objects
                    str(sequence), step=step
                )
                result = self.get_cobs_result(res, False)
                result = (
                    dict(sorted(result.items(), key=lambda x: -x[1])[:limit_number])
                    if limit
                    else result
                )
                result_dict[scheme_path_list[counter]] = result
                first_key, highest_result = next(iter(result.items()))
                highest_results[scheme_path_list[counter]] = {first_key: highest_result}
                counter += 1

        # check if the strain type has sufficient amount of kmer hits
        is_valid = self.has_sufficient_score(highest_results, self.avg_locus_bp_size)
        if not is_valid:
            highest_results["Attention:"] = (
                "This strain type is not reliable due to low kmer hit rates!"
            )
        else:
            handler = PubMLSTHandler()
            # allele_id is of type dict
            flattened = {
                locus: int(list(allele_id.keys())[0].split("_")[-1])
                for locus, allele_id in highest_results.items()
            }
            strain_type_name = handler.get_strain_type_name(flattened, self.scheme_url)
            highest_results["ST_Name"] = strain_type_name
        return [{"Strain type": highest_results}, {"All results": result_dict}]

    def predict(
        self,
        sequence_input: (
            SeqRecord
            | list[SeqRecord]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Path
        ),
        step: int = 1,
        limit: bool = False,
    ) -> MlstResult:
        """
        Get scores for the sequence(s) based on the filters in the model.

        Args:
            sequence_input (Seq): The input sequence for classification
            step (int, optional): The amount of kmers that are passed; defaults to one
            limit (bool, optional): Applying a filter that limits the best result.

        Returns:
            MlstResult: The results of the prediction.

        Raises:
            ValueError: If the sequence input is invalid.
        """
        if isinstance(sequence_input, SeqRecord):
            if sequence_input.id == "<unknown id>":
                sequence_input.id = "test"
            hits = {
                sequence_input.id: self.calculate_hits(sequence_input.seq, step, limit)
            }
            return MlstResult(self.model_display_name, step, hits, None)

        if isinstance(sequence_input, Path):
            return ProbabilisticFilterMlstSchemeModel.predict(
                self,
                get_record_iterator(sequence_input),
                step=step,
                limit=limit,
            )

        if isinstance(
            sequence_input,
            (SeqIO.FastaIO.FastaIterator, SeqIO.QualityIO.FastqPhredIterator),
        ):
            hits = {}
            # individual_seq is a SeqRecord-Object
            for individual_seq in sequence_input:
                individual_hits = self.calculate_hits(individual_seq.seq, step, limit)
                hits[individual_seq.id] = individual_hits
            return MlstResult(self.model_display_name, step, hits, None)
        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a"
            " SeqIO FastaIterator, or a SeqIO FastqPhredIterator"
        )

    def get_cobs_result(
        self,
        cobs_result: cobs_index.SearchResult,
        kmer_threshold: bool,
    ) -> dict:
        """
        Get every entry in a COBS search result.

        Args:
            cobs_result (SearchResult): The result of the prediction.
            kmer_threshold (bool): Applying a kmer threshold to mitigate false positives.

        Returns:
            dict: A dictionary storing the allele id of locus as key and the score as value.
        """
        hits = [
            result for result in cobs_result if not kmer_threshold or result.score > 50
        ]
        return {result.doc_name: result.score for result in hits}

    def sequence_splitter(self, input_sequence: str, allele_len: int) -> list[str]:
        """
        Get an equally divided sequence in form of a list.

        This function is splitting very long sequences into substrings.
        The split is based on sequence and allele length.
        Measures have been taken to not lose kmers while splitting.

        Args:
            input_sequence (str): The sequence of interest.
            allele_len (int): The average length of an allele.

        Returns:
            list[str]: A list containing all substrings of a sequence greater than 10000 bp.

        Raises:
            ValueError: If the sequence input is invalid.
        """

        # An input sequence will have 10000 or more base pairs.
        sequence_len = len(input_sequence)

        if sequence_len < 1000000:
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

    def has_sufficient_score(
        self, highest_results: dict, locus_size: list[int]
    ) -> bool:
        """
        Checks if at least one locus in highest_results has a score >= 0.5 * avg base pair size.

        Args:
            highest_results (dict): Dict where each key is a locus and each value is the kmer score.
            locus_size (list[int]): List of average base pair sizes per locus (in directory order).

        Returns:
            bool: True if any locus score >= 0.5 * its avg base pair size, False otherwise.
        """
        for i, (_, allele_score_dict) in enumerate(highest_results.items()):
            if not allele_score_dict:
                continue  # skip empty values

            # Take the score (the only value) from the nested dict
            score = next(iter(allele_score_dict.values()))
            if score >= 0.5 * locus_size[i]:
                return True
        return False
