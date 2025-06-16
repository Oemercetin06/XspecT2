# Understanding XspecT

## What is XspecT?

XspecT is a tool designed to monitor and characterize pathogens using exact pattern matching of kmers. It allows users to filter for pathogen sequences in metagenomic datasets, classify these sequences on a species level, and perform strain-level typing.

## Key Features
- **Genus-Level Classification**: Classify sequences at the genus level, enabling researchers to quickly identify the presence of specific microbial groups.
- **Species-Level Classification**: Provides detailed classification of sequences at the species level, enhancing the understanding of microbial diversity.
- **Multi-Locus Strain Typing**: Offers the ability to type sequences at the strain level, which is crucial for understanding variations within species.
- **Filtering**: Classification results can be used to filter sequences, enabling analysis of metagenomic samples.
- **Model Management**: XspecT models can be easily downloaded or trained from scratch using the command line interface. Training is possible both from local data, as well as from the NCBI Datasets and PubMLST API.
- **User-friendly Interface**: Next to the command line interface (CLI), a React-based web interface is available for easy interaction and visualization of results.
- **Works with Large Datasets**: Entire folders of input data can be passed to the tool, allowing for efficient processing of large datasets.

## How XspecT Works
At its core, XspecT uses exact pattern matching of kmers to identify and classify sequences. The tool leverages indices of known pathogen sequences stored in XspecT models to match against input data. This process involves:

1. **Kmer Extraction**: The input sequences are processed to extract kmers, which are short sequences of a fixed length.
2. **Pattern Matching**: The extracted kmers are matched against an index of known sequences using exact matching algorithms. The number of matches is recorded, and stored as hits.
3. **Classification**: Based on hits, scores are calculated as the fraction of kmers that match known sequences. These scores are then used to classify the sequences at different taxonomic levels.

### COBS Index
In order to store kmers in a space-efficient manner, XspecT uses a COBS ("Compact Bit-Sliced Signature Index") classic index. This index uses a probabilistic data structure to store kmers, allowing for efficient storage and retrieval. The COBS index is designed to handle large datasets while maintaining fast query performance. More information about the COBS index can be found in the [COBS research paper](https://arxiv.org/abs/1905.09624).