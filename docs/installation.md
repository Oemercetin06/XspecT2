# Installation

## Installation using PyPi

The easiest way to install XspecT is to use PyPi on the latest version of Python:

```bash
pip install xspect
```

Check if the installation was successful:

```bash
xspect --version
```

### Optional Custom Filter Training Requirement

If you would like to train custom filters, you need to install Jellyfish, for example using Conda:

```bash
conda install -c bioconda jellyfish 
```

Note that on Apple Silicon, it is possible that the wrong Jellyfish packages gets installed as binaries are not yet available for ARM. Please refer to the official [Jellyfish project](https://github.com/gmarcais/Jellyfish) for installation guidance.

## Manual Installation

If you would like to manually install XspecT, please clone the Github repository to a local working directory. You can now install XspecT by running:

```bash
pip install .
```

For development purposes, it is recommended to install the package in edit mode:

```bash
pip install -e '.[all]'
```

Check if the installation was successful:

```bash
xspect --version
```

For filter training, please see [Optional Filter Training Requirement]

[Optional Filter Training Requirement]: #optional-custom-filter-training-requirement