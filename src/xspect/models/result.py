"""Module for storing the results of XspecT models."""

from json import dumps
from pathlib import Path


class ModelResult:
    """Class for storing an XspecT model result."""

    def __init__(
        self,
        # we store hits depending on the subsequence as well as on the label
        model_slug: str,
        hits: dict[str, dict[str, int]],
        num_kmers: dict[str, int],
        sparse_sampling_step: int = 1,
        prediction: str = None,
        input_source: str = None,
    ):
        if "total" in hits:
            raise ValueError(
                "'total' is a reserved key and cannot be used as a subsequence"
            )
        self.model_slug = model_slug
        self.hits = hits
        self.num_kmers = num_kmers
        self.sparse_sampling_step = sparse_sampling_step
        self.prediction = prediction
        self.input_source = input_source

    def get_scores(self) -> dict:
        """Return the scores of the model."""
        scores = {
            subsequence: {
                label: round(hits / self.num_kmers[subsequence], 2)
                for label, hits in subsequence_hits.items()
            }
            for subsequence, subsequence_hits in self.hits.items()
        }

        # calculate total scores
        total_num_kmers = sum(self.num_kmers.values())
        total_hits = self.get_total_hits()

        scores["total"] = {
            label: round(hits / total_num_kmers, 2)
            for label, hits in total_hits.items()
        }

        return scores

    def get_total_hits(self) -> dict[str, int]:
        """Return the total hits of the model."""
        total_hits = {label: 0 for label in list(self.hits.values())[0]}
        for _, subseuqence_hits in self.hits.items():
            for label, hits in subseuqence_hits.items():
                total_hits[label] += hits
        return total_hits

    def get_filter_mask(self, label: str, filter_threshold: float) -> dict[str, bool]:
        """Return a mask for filtered subsequences.

        The mask is a dictionary with subsequence names as keys and boolean values
        indicating whether the subsequence is above the filter threshold for the given label.
        A value of -1 for filter_threshold indicates that the subsequence with the maximum score
        for the given label should be returned.
        """
        if filter_threshold < 0 and not filter_threshold == -1 or filter_threshold > 1:
            raise ValueError("The filter threshold must be between 0 and 1.")

        scores = self.get_scores()
        scores.pop("total")
        if not filter_threshold == -1:
            return {
                subsequence: score[label] >= filter_threshold
                for subsequence, score in scores.items()
            }
        else:
            return {
                subsequence: score[label] == max(score.values())
                for subsequence, score in scores.items()
            }

    def get_filtered_subsequence_labels(
        self, label: str, filter_threshold: float = 0.7
    ) -> list[str]:
        """Return the labels of filtered subsequences."""
        return [
            subsequence
            for subsequence, mask in self.get_filter_mask(
                label, filter_threshold
            ).items()
            if mask
        ]

    def to_dict(self) -> dict:
        """Return the result as a dictionary."""
        res = {
            "model_slug": self.model_slug,
            "sparse_sampling_step": self.sparse_sampling_step,
            "hits": self.hits,
            "scores": self.get_scores(),
            "num_kmers": self.num_kmers,
            "input_source": self.input_source,
        }

        if self.prediction is not None:
            res["prediction"] = self.prediction

        return res

    def save(self, path: Path) -> None:
        """Save the result as a JSON file."""
        with open(path, "w", encoding="utf-8") as f:
            f.write(dumps(self.to_dict(), indent=4))
