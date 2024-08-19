"""Project CLI"""

from pathlib import Path
import click
import uvicorn
from xspect import fastapi
from xspect.download_filters import download_test_filters
from xspect.model_management import get_genus_model, get_species_model
from xspect.train import train_ncbi
from xspect.file_io import get_records_by_id


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
@click.argument("path", type=click.Path(exists=True, dir_okay=False, file_okay=True))
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
    """Classify sample(s) from directory PATH."""
    click.echo("Classifying sample...")
    click.echo(f"Step: {step}")
    sequence_input = Path(path)
    species_filter_model = get_species_model(genus)

    if meta:
        genus_filter_model = get_genus_model(genus)
        filter_res = genus_filter_model.predict(sequence_input)
        filtered_sequence_ids = filter_res.get_filtered_subsequences(genus, 0.7)
        print(filtered_sequence_ids)
        sequence_input = get_records_by_id(sequence_input, filtered_sequence_ids)
    res = species_filter_model.predict(sequence_input, step=step)
    print(species_filter_model.display_names[res.prediction])
    print(res.get_scores())


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
