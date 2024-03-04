# How to use the CLI

XspecT comes with a built-in command line interface (CLI), which enables quick classifications without the need to use the web interface. The command line interface can also be used to download and train filters.

After installing XspecT, a list of available commands can be viewed by running:

```bash
xspect --help
```

## Filter downloads

A basic set of pre-trained filters (Acinetobacter and Salonella) can be downloaded using the following command:

```bash
xspect download-filters
```

For the moment, it is not possible to specify exactly which filters should be downloaded.

## Classification

To classify samples, the command

```bash
xspect classify GENUS PATH
```

can be used, when `GENUS` refers to the NCBI genus name of your sample and `PATH` refers to the path to your sample *directory*. This command will classify the species of your sample within the given genus.

The following options are available:

```bash
-s, --species / --no-species    Species classification.
-i, --ic / --no-ic              IC strain typing.
-o, --oxa / --no-oxa            OXA gene family detection.
-m, --meta / --no-meta          Metagenome classification.
-c, --complete                  Use every single k-mer as input for
                                  classification.
-s, --save                      Save results to csv file.
--help                          Show this message and exit.
```

### Species Classification

Species classification is run by default, without the need for further parameters:
```bash
xspect classify Acinetobacter path
```

Species classification can be toggled using the `-s`/`--species` (`--no-species`) option. To run classification without species classification, the option `--no-species` can be used, for example when running a different analysis:

```bash
xspect classify --no-species -i Acinetobacter path
```

### IC Strain Typing

To perform International Clonal (IC) type classification, the `-i`/`--ic` (`--no-ic`) option can be used:

```bash
xspect classify -i Acinetobacter path
```

Please note that IC strain typing is only available for Acinetobacter baumanii.

### OXA Gene Detection

OXA gene detection can be enabled using the `-o`/`--oxa` (`--no-oxa`) option.

```bash
xspect classify -o Acinetobacter path
```

### Metagenome Mode

To analyze a sample in metagenome mode, the `-m`/`--meta` (`--no-meta`) option can be used:

```bash
xspect classify -m Acinetobacter path
```

Compared to normal XspecT modes, this mode first identifies reads belonging to the given genus and continues classification only with the resulting reads and is thus more suitable for metagenomic samples as the resulting runtime is decreased.

## Filter Training

<aside>
⚠️ Depending on genome size and the amount of species, training can take time!

</aside>

In order to train filters, please first ensure [Jellyfish](https://github.com/gmarcais/Jellyfish) is installed.

### NCBI-based filter training

The easiest way to train new filters is to use data from NCBI, which is automatically downloaded and processed by XspecT.

To train a filter with data from NCBI, run the following command:

```bash
xspect train your-ncbi-genus
```

`you-ncbi-genus` can be a genus name from NCBI or an NCBI taxonomy ID.

### Custom data filter training

XspecT filters can also be trained using custom data, which need to be provided as a folder for both filter and SVM training. The provided assembly files need to be in FASTA format and their names should be the species ID and the species name, for example `28901_enterica.fasta`. While the ID can be arbitrary, the standard is NCBI taxon IDs.

The filters can then be trained using:

```bash
xspect train -bf-path directory/1 -svm-path directory/2
```