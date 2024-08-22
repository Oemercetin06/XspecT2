# XspecT - Acinetobacter Species Assignment Tool
![Test](https://github.com/bionf/xspect2/actions/workflows/test.yml/badge.svg)
[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/pylint-dev/pylint)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

<img src="/src/docs/img/logo.png" height="50%" width="50%">
<!-- start intro -->
XspecT is a Python-based tool to taxonomically classify sequence-reads (or assembled genomes) on the species and/or sub-type level using [Bloom Filters] and a [Support Vector Machine]. It also identifies existing [blaOxa-genes] and provides a list of relevant research papers for further information.
<br/><br/>

XspecT utilizes the uniqueness of kmers and compares extracted kmers from the input-data to a reference database. Bloom Filter ensure a fast lookup in this process. For a final prediction the results are classified using a Support Vector Machine. 
<br/>

Local extensions of the reference database are supported.
<br/>

The tool is available as a web-based application and a smaller command line interface.

[Bloom Filters]: https://en.wikipedia.org/wiki/Bloom_filter
[Support Vector Machine]: https://en.wikipedia.org/wiki/Support-vector_machine
[blaOxa-genes]: https://en.wikipedia.org/wiki/Beta-lactamase#OXA_beta-lactamases_(class_D)
<!-- end intro -->

<!-- start quickstart -->
## Installation
To install Xspect, please download the lastest 64 bit Python version and install the package using pip:
```
pip install xspect
```
Please note that Apple Silicon is currently not supported.

## Usage
### Get the Bloomfilters
To download basic pre-trained filters, you can use the built-in command:
```
xspect download-filters
```
Additional species filters can be trained using:
```
xspect train you-ncbi-genus-name
```

### How to run the web app
To run the web app, install and run [XspecT Web](https://github.com/aromberg/xspect-web). Additionally, run XspecT in API mode:
```
xspect api
```

### How to use the XspecT command line interface
Run xspect with the configuration you want to run it with as arguments.
```
xspect classify your-genus path/to/your/input-set
```
For further instructions on how to use the command line interface, please refer to the [documentation] or execute:
```
xspect --help
```
[documentation]: https://bionf.github.io/XspecT2/cli.html
<!-- end quickstart -->