""" Module for storing the results of XspecT models. """


class ModelResult:
    """Class for storing an XspecT model result."""

    def __init__(
        self,
        # we store hits depending on the subsequence as well as on the label
        hits: dict[str, dict[str, int]],
        num_kmers: dict[str, int],
        prediction: str = None,
    ):
        if "total" in hits:
            raise ValueError(
                "'total' is a reserved key and cannot be used as a subsequence"
            )
        self.hits = hits
        self.num_kmers = num_kmers
        self.prediction = prediction

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

    def __dict__(self) -> dict:
        """Return the result as a dictionary."""
        res = {
            "hits": self.hits,
            "scores": self.get_scores(),
            "num_kmers": self.num_kmers,
        }

        if self.prediction is not None:
            res["prediction"] = self.prediction

        return res
