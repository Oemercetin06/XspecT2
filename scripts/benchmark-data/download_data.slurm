#!/bin/bash
#SBATCH --partition=all 
#SBATCH --account=intern              
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8gb             
#SBATCH --job-name="download_acinetobacter"

datasets download genome taxon 469 --filename acinetobacter_dataset.zip --assembly-source GenBank --assembly-version latest --exclude-atypical --dehydrated
unzip -o acinetobacter_dataset.zip -d genomes
datasets rehydrate --directory genomes
rm acinetobacter_dataset.zip

datasets summary taxonomy taxon 469 --rank species --children > aci_species.json