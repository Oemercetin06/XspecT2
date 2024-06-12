"""Project CLI"""

from pathlib import Path
import click
from xspect.download_filters import download_test_filters
from xspect.model_management import get_genus_model, get_species_model
from xspect.train import train_ncbi
from xspect.web_app import app


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
    default=500,
)
def classify(genus, path, meta, step):
    """Classify sample(s) from directory PATH."""
    click.echo("Classifying sample...")
    click.echo(f"Step: {step}")
    sequence_input = Path(path)
    species_filter_model = get_species_model(genus)

    if meta:
        genus_filter_model = get_genus_model(genus)
        filtered_sequences = genus_filter_model.filter(sequence_input)
        prediction, scores = species_filter_model.predict(
            filtered_sequences["Acinetobacter"], step=step
        )
    else:
        prediction, scores = species_filter_model.predict(sequence_input, step=step)

    prediction = species_filter_model.display_names[prediction[0]]

    print(prediction)
    print(scores)


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
    default=500,
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
def web():
    """Open the XspecT web app."""
    port = 8000
    print(f"To open the web app, go to http://localhost:{port}/")
    app.run(host="0.0.0.0", port=port, debug=True, threaded=True)


if __name__ == "__main__":
    cli()
