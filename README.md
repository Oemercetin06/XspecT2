# XspecT
<!-- start intro -->
![Test](https://github.com/bionf/xspect2/actions/workflows/test.yml/badge.svg)
[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/pylint-dev/pylint)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

XspecT is a Python-based tool to taxonomically classify sequence-reads (or assembled genomes) on the species and/or MLST level using [kmer indices] and a [Support Vector Machine].

XspecT utilizes the uniqueness of kmers and compares extracted kmers from the input-data to a kmer index. Probablistic data structures ensure a fast lookup in this process. For a final prediction, the results are classified using a Support Vector Machine.

The tool is available as a web-based application and as a command line interface.

[kmer indices]: https://arxiv.org/abs/1905.09624
[Support Vector Machine]: https://en.wikipedia.org/wiki/Support-vector_machine
<!-- end intro -->

<!-- start quickstart -->
## Installation
To install XspecT, please download the lastest 64 bit Python version and install the package using pip:
```
pip install xspect
```
Please note that Windows and Alpine Linux are currently not supported.

## Usage
### Get the models
To download basic pre-trained models, you can use the built-in command:
```
xspect models download
```
Additional species models can be trained using:
```
xspect models train ncbi
```

### How to run the web app
To run the web app, simply execute:
```
xspect web
```

This will start a local web server. You can access the web app by navigating to `http://localhost:8000` in your web browser.

### How to use the XspecT command line interface
To use the XspecT command line interface, execute `xspect` with the desired subcommand and parameters.

**Example**:
```
xspect classify species
```

If you do not provide the required parameters, the command line interface will prompt you for them.
For further instructions on how to use the command line interface, please refer to the [documentation] or execute:
```
xspect --help
```
[documentation]: https://bionf.github.io/XspecT2/cli/index.html
<!-- end quickstart -->