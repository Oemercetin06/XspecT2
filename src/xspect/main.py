"""Project CLI"""

from pathlib import Path
from uuid import uuid4
from importlib import import_module
import click
from xspect.handlers.pubmlst import PubMLSTHandler
from xspect.model_management import (
    get_available_mlst_schemes,
    get_models,
)

# inline imports lead to "invalid name" issues
# pylint: disable=invalid-name


@click.group()
@click.version_option()
def cli():
    """XspecT CLI."""


@cli.command()
def web():
    """Open the XspecT web application."""
    app = import_module("xspect.web").app
    run = import_module("uvicorn").run

    run(app, host="0.0.0.0", port=8000)


@cli.command(
    name="all",
    help="Run full classification pipeline: genus filtering, species classification, and MLST (if applicable).",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to use.",
    type=click.Choice(get_models().get("Species", [])),
    prompt=True,
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-dir",
    help="Directory for output files (default: auto-generated 'xspect_results_<uuid>' directory).",
    type=click.Path(dir_okay=True, file_okay=False),
    default=None,
)
@click.option(
    "-t",
    "--threshold",
    type=click.FloatRange(0, 1),
    help="Threshold for genus filtering (default: 0.7).",
    default=0.7,
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for '--sparse-sampling-step 500').",
    default=1,
)
@click.option(
    "-n",
    "--display-names",
    help="Includes the display names next to taxonomy-IDs.",
    is_flag=True,
)
@click.option(
    "-v",
    "--validation",
    help="Detects misclassification for small reads or contigs.",
    is_flag=True,
)
def all_pipeline(
    model_genus,
    input_path,
    output_dir,
    threshold,
    sparse_sampling_step,
    display_names,
    validation,
):
    """Run full classification pipeline."""
    import json
    from xspect import filter_sequences, classify
    from xspect.file_io import fasta_endings, fastq_endings

    run_id = uuid4()

    if output_dir is None:
        output_dir = Path(f"xspect_results_{run_id}")
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(exist_ok=True, parents=True)
    input_path = Path(input_path)

    filtered_dir = output_dir / "filtered_sequences"
    filtered_dir.mkdir(exist_ok=True, parents=True)

    genus_filtered_path = filtered_dir / f"genus_filtered_{run_id}.fasta"
    genus_classification_path = output_dir / f"genus_classification_{run_id}.json"
    species_classification_path = output_dir / f"species_classification_{run_id}.json"

    # Step 1: Genus filtering
    click.echo(f"Step 1/3: Filtering for genus {model_genus}...")
    filter_sequences.filter_genus(
        model_genus,
        input_path,
        genus_filtered_path,
        threshold,
        genus_classification_path,
        sparse_sampling_step=sparse_sampling_step,
    )

    ending_wildcards = [f"*.{ending}" for ending in fasta_endings + fastq_endings]
    filtered_files = [p for e in ending_wildcards for p in filtered_dir.glob(e)]

    if not filtered_files:
        click.echo("No sequences passed the genus filter. Pipeline aborted.")
        return

    # Step 2: Species classification on filtered sequences
    click.echo(
        f"Step 2/3: Classifying species for {len(filtered_files)} filtered file(s)..."
    )
    classify.classify_species(
        model_genus,
        filtered_dir,
        species_classification_path,
        sparse_sampling_step,
        display_names,
        validation,
        None,
    )

    species_results = list(output_dir.glob(f"species_classification_{run_id}*.json"))

    # Step 3: Check if we need to run MLST (if any prediction is "470" for abaumannii)
    mlst_needed = False
    for species_result_path in species_results:
        with open(species_result_path, "r", encoding="utf-8") as f:
            species_result = json.load(f)

        prediction = species_result.get("prediction")
        if prediction == "470":
            mlst_needed = True
            click.echo(
                f"Species prediction is 470 (abaumannii) in {species_result_path.name}."
            )

    if mlst_needed:
        click.echo("Step 3/3: Running MLST classification for abaumannii...")

        mlst_schemes = get_available_mlst_schemes()
        if "abaumannii" in mlst_schemes and mlst_schemes["abaumannii"]:
            scheme = mlst_schemes["abaumannii"][0]
            mlst_output_path = output_dir / f"mlst_classification_{run_id}.json"

            classify.classify_mlst(
                filtered_dir,
                "abaumannii",
                scheme,
                mlst_output_path,
                False,
            )
            click.echo(f"MLST classification completed: {mlst_output_path.name}")
        else:
            click.echo(
                "Warning: No MLST schemes available for abaumannii. Skipping MLST classification."
            )
    else:
        click.echo(
            "Step 3/3: Not running MLST classification (organism is not Acinetobacter baumannii)."
        )

    click.echo("\nPipeline completed successfully!")
    click.echo(f"Results saved in: {output_dir}")


# # # # # # # # # # # # # # #
# Model management commands #
# # # # # # # # # # # # # # #
@cli.group()
def models():
    """Model management commands."""


@models.command(
    help="Download models from the internet.",
)
def download():
    """Download models."""
    click.echo("Downloading models, this may take a while...")
    download_test_models = import_module("xspect.download_models").download_test_models
    download_test_models()


@models.command(
    name="list",
    help="List all models in the model directory.",
)
def list_models():
    """List models."""
    available_models = get_models()
    if not available_models:
        click.echo("No models found.")
        return
    click.echo("Models found:")
    click.echo("--------------")
    for model_type, names in available_models.items():
        if not names:
            continue
        click.echo(f"  {model_type}:")
        for name in names:
            click.echo(f"    - {name}")


@models.group()
def train():
    """Train models."""


@train.command(
    name="ncbi",
    help="Train a species and a genus model based on NCBI data.",
)
@click.option("-g", "--genus", "model_genus", prompt=True)
@click.option("--svm_steps", type=int, default=1)
@click.option(
    "--author",
    help="Author of the model.",
    default=None,
)
@click.option(
    "--author-email",
    help="Email of the author.",
    default=None,
)
@click.option(
    "--min-n50",
    type=int,
    help="Minimum contig N50 to filter the accessions (default: 10000).",
    default=10000,
)
@click.option(
    "--include-atypical/--exclude-atypical",
    help="Include or exclude atypical accessions (default: exclude).",
    default=False,
)
@click.option(
    "--allow-inconclusive",
    is_flag=True,
    help="Allow the use of accessions with inconclusive taxonomy check status for training.",
    default=False,
)
@click.option(
    "--allow-candidatus",
    is_flag=True,
    help="Allow the use of Candidatus species for training.",
    default=False,
)
@click.option(
    "--allow-sp",
    is_flag=True,
    help="Allow the use of species with 'sp.' in their names for training.",
    default=False,
)
def train_ncbi(
    model_genus,
    svm_steps,
    author,
    author_email,
    min_n50,
    include_atypical,
    allow_inconclusive,
    allow_candidatus,
    allow_sp,
):
    """Train a species and a genus model based on NCBI data."""
    click.echo(f"Training {model_genus} species and genus metagenome model.")
    try:
        train_from_ncbi = import_module("xspect.train").train_from_ncbi

        train_from_ncbi(
            model_genus,
            svm_steps,
            author,
            author_email,
            min_n50=min_n50,
            exclude_atypical=not include_atypical,
            allow_inconclusive=allow_inconclusive,
            allow_candidatus=allow_candidatus,
            allow_sp=allow_sp,
        )
    except ValueError as e:
        click.echo(f"Error: {e}")
        return
    click.echo(f"Training of {model_genus} model finished.")


@train.command(
    name="directory",
    help="Train a species (and possibly a genus) model based on local data.",
)
@click.option("-g", "--genus", "model_genus", prompt=True)
@click.option(
    "-i",
    "--input-path",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
)
@click.option(
    "--meta",
    is_flag=True,
    help="Train a metagenome model for the genus.",
    default=True,
)
@click.option(
    "--svm-steps",
    type=int,
    help="SVM Sparse sampling step size (e. g. only every 500th kmer for step=500).",
    default=1,
)
@click.option(
    "--author",
    help="Author of the model.",
    default=None,
)
@click.option(
    "--author-email",
    help="Email of the author.",
    default=None,
)
def train_directory(model_genus, input_path, svm_steps, meta, author, author_email):
    """Train a model based on data from a directory for a given genus."""
    click.echo(f"Training {model_genus} model with {svm_steps} SVM steps.")
    train_from_directory = import_module("xspect.train").train_from_directory

    train_from_directory(
        model_genus,
        Path(input_path),
        svm_step=svm_steps,
        meta=meta,
        author=author,
        author_email=author_email,
    )


@train.command(
    name="mlst",
    help="Train a MLST model based on PubMLST data.",
)
@click.option(
    "--organism", "organism", help="Underlying organism for the MLST model.", type=str
)
@click.option(
    "--mlst-scheme",
    "scheme",
    help="MLST scheme to use for the model.",
    type=str,
)
@click.option(
    "--author",
    help="Author of the model.",
    default=None,
)
@click.option(
    "--author-email",
    help="Email of the author.",
    default=None,
)
def train_mlst(organism, scheme, author, author_email):
    """Download alleles and train MLST models."""
    handler = PubMLSTHandler()
    available_organisms = handler.get_available_organisms()
    if not organism:
        organism = click.prompt(
            "Please enter the organism you want to train the MLST model for:",
            type=click.Choice(available_organisms),
        )
    elif organism not in available_organisms:
        raise click.BadParameter(
            f"Organism '{organism}' not found. Available organisms: {', '.join(available_organisms)}"
        )

    available_schemas = handler.get_available_schemes(organism)
    if scheme:
        if scheme not in available_schemas:
            raise click.BadParameter(
                f"Scheme '{scheme}' not found for organism '{organism}'. "
                f"Available schemes: {', '.join(available_schemas)}"
            )
    else:
        scheme = click.prompt(
            "Please enter the scheme you want to train the MLST model for:",
            type=click.Choice(available_schemas),
        )

    train_mlst_model = import_module("xspect.train").train_mlst
    train_mlst_model(organism, scheme, author, author_email)


# # # # # # # # # # # # # # #
# Classification commands   #
# # # # # # # # # # # # # # #
@cli.group(
    name="classify",
    help="Classify sequences using XspecT models.",
)
def classify_seqs():
    """Classification commands."""


@classify_seqs.command(
    name="genus",
    help="Classify samples using a genus model.",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to classify.",
    type=click.Choice(get_models().get("Genus", [])),
    prompt=True,
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    default=Path(".") / f"result_{uuid4()}.json",
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for '--sparse-sampling-step 500').",
    default=1,
)
def classify_genus(model_genus, input_path, output_path, sparse_sampling_step):
    """Classify samples using a genus model."""
    click.echo("Classifying...")
    classify = import_module("xspect.classify")

    classify.classify_genus(
        model_genus, Path(input_path), Path(output_path), sparse_sampling_step
    )


@classify_seqs.command(
    name="species",
    help="Classify samples using a species model.",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to classify.",
    type=click.Choice(get_models().get("Species", [])),
    prompt=True,
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    default=Path(".") / f"result_{uuid4()}.json",
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for '--sparse-sampling-step 500').",
    default=1,
)
@click.option(
    "-n",
    "--display-names",
    help="Includes the display names next to taxonomy-IDs.",
    is_flag=True,
)
@click.option(
    "-v",
    "--validation",
    help="Detects misclassification for small reads or contigs.",
    is_flag=True,
)
@click.option(
    "--exclude-species",
    help="Comma-separated list of species IDs to exclude from classification.",
    type=str,
    default=None,
)
def classify_species(
    model_genus,
    input_path,
    output_path,
    sparse_sampling_step,
    display_names,
    validation,
    exclude_species,
):
    """Classify samples using a species model."""
    click.echo("Classifying...")
    classify = import_module("xspect.classify")

    exclude_ids = None
    if exclude_species:
        exclude_ids = [species_id.strip() for species_id in exclude_species.split(",")]

    classify.classify_species(
        model_genus,
        Path(input_path),
        Path(output_path),
        sparse_sampling_step,
        display_names,
        validation,
        exclude_ids,
    )


@classify_seqs.command(
    name="mlst",
    help="Classify samples using a MLST model.",
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA-file for mlst identification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "--organism",
    "organism",
    help="Underlying organism for the MLST model.",
    type=click.Choice(get_available_mlst_schemes().keys()),
    prompt=True,
)
@click.option(
    "--mlst-scheme",
    "mlst_scheme",
    help="MLST scheme to use.",
    type=str,
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    default=Path(".") / f"MLST_result_{uuid4()}.json",
)
@click.option(
    "-l", "--limit", is_flag=True, help="Limit the output to 5 results for each locus."
)
def classify_mlst(input_path, organism, mlst_scheme, output_path, limit):
    """MLST classify a sample."""
    mlst_schemes = get_available_mlst_schemes()
    if not mlst_scheme:
        mlst_scheme = click.prompt(
            "Please enter the MLST scheme you want to use:",
            type=click.Choice(mlst_schemes[organism]),
        )
    elif mlst_scheme not in mlst_schemes.get(organism, []):
        raise click.BadParameter(
            f"Scheme '{mlst_scheme}' not found for organism '{organism}'. "
            f"Available schemes: {', '.join(mlst_schemes.get(organism, []))}"
        )

    click.echo("Classifying...")
    classify = import_module("xspect.classify")
    classify.classify_mlst(
        Path(input_path), organism, mlst_scheme, Path(output_path), limit
    )


# # # # # # # # # # # # # # #
# Filtering commands        #
# # # # # # # # # # # # # # #
@cli.group(
    name="filter",
    help="Filter sequences using XspecT models.",
)
def filter_seqs():
    """Filter commands."""


@filter_seqs.command(
    name="genus",
    help="Filter sequences using a genus model.",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to use for filtering.",
    type=click.Choice(get_models().get("Species", [])),
    prompt=True,
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    prompt=True,
    default=Path(".") / f"genus_filtered_{uuid4()}.fasta",
)
@click.option(
    "--classification-output-path",
    help="Optional path to the classification output file.",
    type=click.Path(dir_okay=False, file_okay=True),
)
@click.option(
    "-t",
    "--threshold",
    type=click.FloatRange(0, 1),
    help="Threshold for filtering (default: 0.7).",
    default=0.7,
    prompt=True,
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for '--sparse-sampling-step 500').",
    default=1,
)
def filter_genus(
    model_genus,
    input_path,
    output_path,
    classification_output_path,
    threshold,
    sparse_sampling_step,
):
    """Filter samples using a genus model."""
    click.echo("Filtering...")
    filter_sequences = import_module("xspect.filter_sequences")

    filter_sequences.filter_genus(
        model_genus,
        Path(input_path),
        Path(output_path),
        threshold,
        Path(classification_output_path) if classification_output_path else None,
        sparse_sampling_step=sparse_sampling_step,
    )


@filter_seqs.command(
    name="species",
    help="Filter sequences using a species model.",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to use for filtering.",
    type=click.Choice(get_models().get("Species", [])),
    prompt=True,
)
@click.option(
    "-s",
    "--species",
    "model_species",
    help="Species of the model to filter for.",
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    prompt=True,
    default=Path(".") / f"species_filtered_{uuid4()}.fasta",
)
@click.option(
    "--classification-output-path",
    help="Optional path to the classification output file.",
    type=click.Path(dir_okay=False, file_okay=True),
)
@click.option(
    "-t",
    "--threshold",
    type=float,
    help="Threshold for filtering (default: 0.7). Use -1 to filter for the highest scoring "
    "species.",
    default=0.7,
    prompt=True,
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for "
    "'--sparse-sampling-step 500').",
    default=1,
)
def filter_species(
    model_genus,
    model_species,
    input_path,
    output_path,
    threshold,
    classification_output_path,
    sparse_sampling_step,
):
    """Filter a sample using the species model."""

    if threshold != -1 and (threshold < 0 or threshold > 1):
        raise click.BadParameter(
            "Threshold must be between 0 and 1, or -1 for filtering by the highest "
            "scoring species."
        )

    get_model_metadata = import_module("xspect.model_management").get_model_metadata

    available_species = get_model_metadata(f"{model_genus}-species")["display_names"]
    available_species = {
        id: name.replace(f"{model_genus} ", "")
        for id, name in available_species.items()
    }
    if not model_species:
        sorted_available_species = sorted(available_species.values())
        model_species = click.prompt(
            f"Please enter the species name: {model_genus}",
            type=click.Choice(sorted_available_species, case_sensitive=False),
        )
    if model_species not in available_species.values():
        raise click.BadParameter(
            f"Species '{model_species}' not found in the {model_genus} species model."
        )

    # get the species ID from the name
    model_species = [
        id
        for id, name in available_species.items()
        if name.lower() == model_species.lower()
    ][0]

    click.echo("Filtering...")
    filter_sequences = import_module("xspect.filter_sequences")

    filter_sequences.filter_species(
        model_genus,
        model_species,
        Path(input_path),
        Path(output_path),
        threshold,
        Path(classification_output_path) if classification_output_path else None,
        sparse_sampling_step=sparse_sampling_step,
    )


if __name__ == "__main__":
    cli()
