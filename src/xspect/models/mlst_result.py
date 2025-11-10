""" "Module for storing MLST results."""

import json
from pathlib import Path


class MlstResult:
    """Class for storing MLST results."""

    def __init__(
        self,
        scheme_model: str,
        steps: int,
        hits: dict[str, list[dict]],
        input_source: str | None = None,
    ):
        """Initialise an MlstResult object."""
        self.scheme_model = scheme_model
        self.steps = steps
        self.hits = hits
        self.input_source = input_source

    def get_results(self) -> dict:
        """
        Stores the result of a prediction in a dictionary.

        Returns:
            dict: The result dictionary with s sequence ID as key and the Strain type as value.
        """
        return dict(self.hits.items())

    def to_dict(self) -> dict:
        """
        Converts all attributes into one dictionary.

        Returns:
            dict: The dictionary containing all metadata of a run.
        """
        result = {
            "Scheme": self.scheme_model,
            "Steps": self.steps,
            "Results": self.get_results(),
            "Input_source": self.input_source,
        }
        return result

    def save(self, output_path: Path | str) -> None:
        """
        Saves the result as a JSON file.

        Args:
            output_path (Path,str): The path where the results are saved.
        """

        if isinstance(output_path, str):
            output_path = Path(output_path)

        output_path.parent.mkdir(exist_ok=True, parents=True)
        json_object = json.dumps(self.to_dict(), indent=4)

        with open(output_path, "w", encoding="utf-8") as file:
            file.write(json_object)
