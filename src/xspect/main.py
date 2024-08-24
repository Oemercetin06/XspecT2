"""Project CLI"""

from pathlib import Path
import datetime
import click
import uvicorn
from xspect import fastapi
from xspect.download_filters import download_test_filters
from xspect.train import train_ncbi
from xspect.models.result import (
    StepType,
)
from xspect.definitions import get_xspect_runs_path, fasta_endings, fastq_endings
from xspect.pipeline import ModelExecution, Pipeline, PipelineStep


@click.group()
@click.version_option()
def cli():
    """XspecT CLI."""


@cli.command()
def download_filters():
    """Download filters."""
    click.echo("Downloading filters, this may take a while...")
    download_test_filters("https://xspect2.s3.eu-central-1.amazonaws.com/models.zip")


@cli.command()
@click.argument("genus")
@click.argument("path", type=click.Path(exists=True, dir_okay=True, file_okay=True))
@click.option(
    "-m",
    "--meta/--no-meta",
    help="Metagenome classification.",
    default=False,
)
@click.option(
    "-s",
    "--step",
    help="Sparse sampling step size (e. g. only every 500th kmer for step=500).",
    default=1,
)
def classify(genus, path, meta, step):
    """Classify sample(s) from file or directory PATH."""
    click.echo("Classifying...")
    click.echo(f"Step: {step}")

    file_paths = []
    if Path(path).is_dir():
        file_paths = [
            f
            for f in Path(path).iterdir()
            if f.is_file() and f.suffix[1:] in fasta_endings + fastq_endings
        ]
    else:
        file_paths = [Path(path)]

    # define pipeline
    pipeline = Pipeline(genus + " classification", "Test Author", "test@example.com")
    species_execution = ModelExecution(
        genus.lower() + "-species", sparse_sampling_step=step
    )
    if meta:
        species_filtering_step = PipelineStep(
            StepType.FILTERING, genus, 0.7, species_execution
        )
        genus_execution = ModelExecution(
            genus.lower() + "-genus", sparse_sampling_step=step
        )
        genus_execution.add_pipeline_step(species_filtering_step)
        pipeline.add_pipeline_step(genus_execution)
    else:
        pipeline.add_pipeline_step(species_execution)

    for idx, file_path in enumerate(file_paths):
        run = pipeline.run(file_path)
        time_str = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        save_path = get_xspect_runs_path() / f"run_{time_str}.json"
        run.save(save_path)
        print(
            f"[{idx+1}/{len(file_paths)}] Run finished. Results saved to '{save_path}'."
        )


@cli.command()
@click.argument("genus")
@click.option(
    "-bf-path",
    "--bf-assembly-path",
    help="Path to assembly directory for Bloom filter training.",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
@click.option(
    "-svm-path",
    "--svm-assembly-path",
    help="Path to assembly directory for SVM training.",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
@click.option(
    "-s",
    "--svm-step",
    help="SVM Sparse sampling step size (e. g. only every 500th kmer for step=500).",
    default=1,
)
def train(genus, bf_assembly_path, svm_assembly_path, svm_step):
    """Train model."""

    if bf_assembly_path or svm_assembly_path:
        raise NotImplementedError(
            "Training with specific assembly paths is not yet implemented."
        )
    try:
        train_ncbi(genus, svm_step=svm_step)
    except ValueError as e:
        raise click.ClickException(str(e)) from e


@cli.command()
def api():
    """Open the XspecT FastAPI."""
    uvicorn.run(fastapi.app, host="0.0.0.0", port=8000)


if __name__ == "__main__":
    cli()
