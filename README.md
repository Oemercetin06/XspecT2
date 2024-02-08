# XspecT - Acinetobacter Species Assignment Tool
<img src="/src/xspect/static/Logo.png" height="50%" width="50%">
<!-- start intro -->
XspecT is a Python-based tool to taxonomically classify sequence-reads (or assembled genomes) on the species and/or sub-type level using [Bloom Filters](https://en.wikipedia.org/wiki/Bloom_filter) and a [Support Vector Machine](https://en.wikipedia.org/wiki/Support-vector_machine). It also identifies existing [blaOxa-genes](https://en.wikipedia.org/wiki/Beta-lactamase#OXA_beta-lactamases_(class_D)) and provides a list of relevant research papers for further information.
<br/><br/>

XspecT utilizes the uniqueness of kmers and compares extracted kmers from the input-data to a reference database. Bloom Filter ensure a fast lookup in this process. For a final prediction the results are classified using a Support Vector Machine. 
<br/>

Local extensions of the reference database are supported.
<br/>

The tool is available as a web-based application and a smaller command line interface.
<!-- end intro -->

<!-- start quickstart -->
## Installation
To install Xspect, please download the lastest 64 bit Python version and install the package using pip:
```
pip install xspect
```
If you would like to train filters yourself, you need to install Jellyfish, which is used to count distinct k-meres in the assemblies. It can be installed using bioconda:
```
conda install -c bioconda jellyfish
```
On Apple Silicon, it is possible that this command installs an incorrect Jellyfish package. Please refer to the official [Jellyfish project](https://github.com/gmarcais/Jellyfish) for installation guidance.

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
Run the following command lines in a console, a browser window will open automatically after the application is fully loaded.
```
xspect web
```

### How to use the XspecT command line interface
Run xspect with the configuration you want to run it with as arguments.
```
xspect classify your-genus path/to/your/input-set
```
For further instructions on how to use the command line interface, execute:
```
xspect --help
```
<!-- end quickstart -->

## Input Data
XspecT is able to use either raw sequence-reads (FASTQ-format .fq/.fastq) or already assembled genomes (FASTA-format .fasta/.fna). Using sequence-reads saves up the assembly process but high-quality reads with a low error-rate are needed (e.g. Illumina-reads).

The amount of reads that will be used has to be set by the user when using sequence-reads. The minimum amount is 5000 reads for species classification and 500 reads for sub-type classification. The maximum number of reads is limited by the browser and is usually around ~8 million reads. Using more reads will lead to a increased runtime (xsec./1mio reads).
