""" Module with XspecT global run class, which summarizes individual model results. """

import json
from pathlib import Path
from xspect.models.result import ModelResult


class Run:
    """Class for storing the results of an XspecT run."""

    def __init__(self, display_name: str, input_file: str):
        self.display_name = display_name
        self.input_file = input_file
        self.results = []

    def add_result(self, result: ModelResult):
        """Add a result to the run."""
        self.results.append(result)

    def to_dict(self) -> dict:
        """Return the run as a dictionary."""
        return {
            "display_name": self.display_name,
            "input_file": str(self.input_file),
            "results": (
                [result.to_dict() for result in self.results] if self.results else []
            ),
        }

    def to_json(self) -> str:
        """Return the run as a JSON string."""
        json_dict = self.to_dict()
        return json.dumps(json_dict, indent=4)

    def save(self, path: Path) -> None:
        """Save the run as a JSON file."""
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.to_json())
