"""Project CLI"""

import webbrowser
import click
from xspect.XspecT_mini import xspecT_mini
from xspect.download_filters import download_test_filters
from xspect.XspecT_trainer import train as x_train
from xspect.WebApp import app


@click.group()
@click.version_option()
def cli():
    """XspecT CLI."""


@cli.command()
def download_filters():
    """Download filters."""
    click.echo("Downloading filters, this may take a while...")
    download_test_filters(
        "https://applbio.biologie.uni-frankfurt.de/download/xspect/filters.zip"
    )


# todo: add read amount option -> why 342480?
@cli.command()
@click.argument("genus")
@click.argument("path", type=click.Path(exists=True, dir_okay=True, file_okay=False))
@click.option(
    "-s", "--species/--no-species", help="Species classification.", default=True
)
@click.option("-i", "--ic/--no-ic", help="IC strain typing.", default=False)
@click.option("-o", "--oxa/--no-oxa", help="OXA gene family detection.", default=False)
@click.option(
    "-m",
    "--metagenome/--no-metagenome",
    help="Metagenome classification.",
    default=False,
)
@click.option(
    "-c",
    "--complete",
    help="Use every single k-mer as input for classification.",
    is_flag=True,
    default=False,
)
@click.option(
    "-s", "--save", help="Save results to csv file.", is_flag=True, default=False
)
def classify(genus, path, species, ic, oxa, metagenome, complete, save):
    """Classify sample(s) from directory PATH."""
    click.echo("Classifying sample...")
    mode = 500
    if complete:
        mode = 1
    file_format = "fasta"
    read_amount = 342480

    xspecT_mini(
        path,
        species,
        ic,
        oxa,
        file_format,
        read_amount,
        save,
        metagenome,
        genus,
        mode,
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
    "-c",
    "--complete",
    help="Train filter on every single k-mer.",
    is_flag=True,
    default=False,
)
@click.option(
    "--check",
    help="Check if metagenome file was correctly created.",
    is_flag=True,
    default=False,
)
def train(genus, bf_assembly_path, svm_assembly_path, complete, check):
    """Train model."""
    mode = "1"
    if bf_assembly_path and svm_assembly_path:
        mode = "2"
    if check:
        mode = "3"
    x_train(genus, mode, complete, bf_assembly_path, svm_assembly_path, "")


@cli.command()
def web():
    """Open the XspecT web app."""
    webbrowser.open("http://localhost:8000")
    app.run(host="0.0.0.0", port=8000, debug=True, threaded=True)


if __name__ == "__main__":
    cli()
