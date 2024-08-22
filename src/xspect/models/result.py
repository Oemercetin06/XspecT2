""" Module for storing the results of XspecT models. """

from enum import Enum


def get_last_processing_step(result: "ModelResult") -> "ModelResult":
    """Get the last subprocessing step of the result. First path only."""
    last_step = result
    while last_step.subprocessing_steps:
        last_step = last_step.subprocessing_steps[-1].result
    return last_step


class StepType(Enum):
    """Enum for defining the type of a subprocessing step."""

    PREDICTION = 1
    FILTERING = 2

    def __str__(self) -> str:
        return self.name.lower()


class SubprocessingStep:
    """Class for storing a subprocessing step of an XspecT model."""

    def __init__(
        self,
        subprocessing_type: StepType,
        label: str,
        treshold: float,
        result: "ModelResult",
    ):
        self.subprocessing_type = subprocessing_type
        self.label = label
        self.treshold = treshold
        self.result = result

    def to_dict(self) -> dict:
        """Return the subprocessing step as a dictionary."""
        return {
            "subprocessing_type": str(self.subprocessing_type),
            "label": self.label,
            "treshold": self.treshold,
            "result": self.result.to_dict() if self.result else {},
        }


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
        self.subprocessing_steps = []

    def add_subprocessing_step(self, subprocessing_step: SubprocessingStep) -> None:
        """Add a subprocessing step to the result."""
        if subprocessing_step.label in self.subprocessing_steps:
            raise ValueError(
                f"Subprocessing step {subprocessing_step.label} already exists in the result"
            )
        self.subprocessing_steps.append(subprocessing_step)

    def get_scores(self) -> dict:
        """Return the scores of the model."""
        scores = {
            subsequence: {
                label: round(hits / self.num_kmers[subsequence], 2)
                for label, hits in subseuqence_hits.items()
            }
            for subsequence, subseuqence_hits in self.hits.items()
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
        """Return a mask for filtered subsequences."""
        if filter_threshold < 0 or filter_threshold > 1:
            raise ValueError("The filter threshold must be between 0 and 1.")

        scores = self.get_scores()
        scores.pop("total")
        return {
            subsequence: score[label] >= filter_threshold
            for subsequence, score in scores.items()
        }

    def get_filtered_subsequences(self, label: str, filter_threshold: 0.7) -> list[str]:
        """Return the filtered subsequences."""
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
            "subprocessing_steps": [
                subprocessing_step.to_dict()
                for subprocessing_step in self.subprocessing_steps
            ],
        }

        if self.prediction is not None:
            res["prediction"] = self.prediction

        return res
