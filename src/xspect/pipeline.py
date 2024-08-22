""" Module for defining the Pipeline class. """

import json
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from xspect.file_io import get_records_by_id
from xspect.models.result import StepType, SubprocessingStep
from xspect.run import Run
from xspect.models.result import ModelResult
from xspect.model_management import get_model_by_slug


class ModelExecution:
    """Class for storing a processing step of an XspecT pipeline."""

    def __init__(
        self,
        model_slug: str,
        sparse_sampling_step: int = 1,
    ):
        self.model_slug = model_slug
        self.sparse_sampling_step = sparse_sampling_step
        self.pipeline_steps = []

    def add_pipeline_step(
        self,
        pipeline_step: "PipelineStep",
    ):
        """Add a subprocessing step to the pipeline step."""
        self.pipeline_steps.append(pipeline_step)

    def to_dict(self) -> dict:
        """Return the processing step as a dictionary."""
        return {
            "model_slug": self.model_slug,
            "sparse_sampling_step": self.sparse_sampling_step,
            "pipeline_steps": [
                pipeline_step.to_dict() for pipeline_step in self.pipeline_steps
            ],
        }

    def run(
        self,
        sequence_input: (
            SeqRecord
            | list[SeqRecord]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Path
        ),
    ) -> ModelResult:
        """Run the model on a given input."""
        model = get_model_by_slug(self.model_slug)
        model_result = model.predict(sequence_input, step=self.sparse_sampling_step)

        for pipeline_step in self.pipeline_steps:
            if pipeline_step.subprocessing_type == StepType.PREDICTION:
                score = model_result.get_scores()["total"][pipeline_step.label]
                if score >= pipeline_step.treshold:
                    prediction_model_result = pipeline_step.model_execution.run(
                        sequence_input
                    )
                    subprocessing_step = SubprocessingStep(
                        pipeline_step.subprocessing_type,
                        pipeline_step.label,
                        pipeline_step.treshold,
                        prediction_model_result,
                    )
                    model_result.add_subprocessing_step(subprocessing_step)
            elif pipeline_step.subprocessing_type == StepType.FILTERING:
                filtered_sequence_ids = model_result.get_filtered_subsequences(
                    pipeline_step.label, pipeline_step.treshold
                )
                sequence_input = get_records_by_id(
                    sequence_input, filtered_sequence_ids
                )

                filtering_model_result = None
                if sequence_input:
                    filtering_model_result = pipeline_step.model_execution.run(
                        sequence_input
                    )

                subprocessing_step = SubprocessingStep(
                    pipeline_step.subprocessing_type,
                    pipeline_step.label,
                    pipeline_step.treshold,
                    filtering_model_result,
                )
                model_result.add_subprocessing_step(subprocessing_step)
            else:
                raise ValueError(
                    f"Invalid subprocessing type {pipeline_step.subprocessing_type}"
                )

        return model_result


class PipelineStep:
    """Class for storing a subprocessing step of an XspecT model."""

    def __init__(
        self,
        subprocessing_type: StepType,
        label: str,
        treshold: float,
        model_execution: ModelExecution,
    ):
        self.subprocessing_type = subprocessing_type
        self.label = label
        self.treshold = treshold
        self.model_execution = model_execution

    def to_dict(self) -> dict:
        """Return the subprocessing step as a dictionary."""
        return {
            "subprocessing_type": str(self.subprocessing_type),
            "label": self.label,
            "treshold": self.treshold,
            "model_execution": self.model_execution.to_dict(),
        }


class Pipeline:
    """Class for storing an XspecT pipeline consisting of multiple model processing steps."""

    def __init__(self, display_name: str, author: str, author_email: str):
        self.display_name = display_name
        self.author = author
        self.author_email = author_email
        self.model_executions = []

    def add_pipeline_step(
        self,
        pipeline_step: ModelExecution,
    ):
        """Add a processing step to the pipeline."""
        self.model_executions.append(pipeline_step)

    def to_dict(self) -> dict:
        """Return the pipeline as a dictionary."""
        return {
            "display_name": self.display_name,
            "author": self.author,
            "author_email": self.author_email,
            "model_executions": [
                model_execution.to_dict() for model_execution in self.model_executions
            ],
        }

    def to_json(self) -> str:
        """Return the pipeline as a JSON string."""
        return json.dumps(self.to_dict())

    def save(self, path: Path) -> None:
        """Save the pipeline as a JSON file."""
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.to_json())

    @staticmethod
    def from_file(path: Path) -> "Pipeline":
        """Load the pipeline from a JSON file."""
        with open(path, "r", encoding="utf-8") as f:
            pipeline_json = json.load(f)
            pipeline = Pipeline(
                pipeline_json["display_name"],
                pipeline_json["author"],
                pipeline_json["author_email"],
            )
            for model_execution in pipeline_json["model_executions"]:
                model_execution = ModelExecution(
                    model_execution["model_slug"],
                    model_execution["sparse_sampling_step"],
                )
                for pipeline_step in model_execution["pipeline_steps"]:
                    model_execution.add_pipeline_step(
                        PipelineStep(
                            StepType(pipeline_step["subprocessing_type"]),
                            pipeline_step["label"],
                            pipeline_step["treshold"],
                            ModelExecution(
                                pipeline_step["model_execution"]["model_slug"],
                                pipeline_step["model_execution"][
                                    "sparse_sampling_step"
                                ],
                            ),
                        )
                    )
                pipeline.add_pipeline_step(model_execution)
            return pipeline

    def run(self, input_file: Path) -> Run:
        """Run the pipeline on a given input."""
        run = Run(self.display_name, input_file)

        for model_execution in self.model_executions:
            result = model_execution.run(input_file)
            run.add_result(result)

        return run
