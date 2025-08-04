# Benchmark

To benchmark XspecT performance, you can use the Nextflow workflow provided in the `scripts/benchmark` directory. This workflow allows you to run XspecT on a set of samples and measure species classification accuracy on both whole genomes, as well as on simulated reads.

Before you run the benchmark, you first need to get benchmarking data, for example from NCBI. To do so, you can use the bash script in `scripts/benchmark-data` to download the data using the [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/), which needs to be installed first. The script will download all available Acinetobacter genomes, as well as taxonomic data.

To run the benchmark, install [Nextflow](https://www.nextflow.io/docs/latest/install.html) and run the following command:

```bash
nextflow run scripts/benchmark
```

This will execute the benchmark workflow, which will classify the samples, as well as reads generated from them, using XspecT. The results will be saved in the `results` directory:

- `results/classifications.tsv` for the classifications of the whole genomes
- `results/read_classifications.tsv` for the classifications of the simulated reads
- `results/confusion_matrix.png` for the confusion matrix of whole genome classifications
- `results/mismatches_confusion_matrix.png` for a confusion matrix filtered on mismatches of whole genome classifications
- `results/stats.txt` for the statistics of the benchmark run