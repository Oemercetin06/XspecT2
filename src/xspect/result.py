""" Module with XspecT global result class, which summarizes individual model results. """

import xspect.models.model_result


class Result:
    """Class for storing an XspecT result."""

    def __init__(self, display_name: str, input_file: str):
        self.display_name = display_name
        self.input_file = input_file
        self.processing_steps = {}

    def add_processing_step(self, step_name: str, results: dict[str, ModelResult]):
        """Add a processing step to the result."""
        if step_name in self.processing_steps:
            raise ValueError(f"Step {step_name} already exists in the result")
        self.processing_steps[step_name] = {"result": results}

    def add_subprocessing_step(
        self,
        preprocessing_steps: list[str],
        substep_name: str,
        results: dict[str, ModelResult],
    ):
        """Add a subprocessing step to the result."""
        try:
            last_step = self.processing_steps[preprocessing_steps[0]]
            for step in preprocessing_steps[1:]:
                last_step = last_step[step]
        except KeyError as e:
            raise ValueError(
                f"Step {preprocessing_steps} does not exist in the result"
            ) from e
        if substep_name in last_step:
            raise ValueError(f"Substep {substep_name} already exists in the result")
        last_step[substep_name] = {"result": results}
