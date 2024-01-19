"""Project CLI"""

import click

from download_filters import download_test_filters


@click.group
def cli():
    pass


@cli.command()
def download_filters():
    """Download filters."""
    click.echo("Downloading filters...")
    download_test_filters(
        "https://xspect.s3.eu-central-1.amazonaws.com/test_filters.zip"
    )


if __name__ == "__main__":
    cli()
