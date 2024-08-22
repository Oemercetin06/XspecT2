from pathlib import Path
from xspect.pipeline import Pipeline, PipelineStep, ModelExecution
from xspect.models.result import StepType


def test_save_and_load(tmpdir):
    """Test saving and loading a pipeline."""
    path = Path(tmpdir) / "test_pipeline.json"
    pipeline = Pipeline("Sample Pipeline", "John Doe", "john.doe@example.com")

    model_execution2 = ModelExecution("test_model2")
    pipeline_step = PipelineStep(
        StepType.FILTERING, "test_label", 0.7, model_execution2
    )
    model_execution = ModelExecution("test_model")
    model_execution.add_pipeline_step(pipeline_step)

    pipeline.save(path)

    loaded_pipeline = Pipeline.from_file(path)

    assert loaded_pipeline.display_name == pipeline.display_name
    assert loaded_pipeline.author == pipeline.author
    assert loaded_pipeline.author_email == pipeline.author_email
    assert len(loaded_pipeline.model_executions) == len(pipeline.model_executions)
    assert loaded_pipeline.to_dict() == pipeline.to_dict()
