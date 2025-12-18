#!/usr/bin/env nextflow

include { classifySample as classifyAssembly } from './classify'
include { classifySample as classifyRead } from './classify'
include { strain_species_mapping } from '../nextflow-utils'
include { confusionMatrix as assemblyConfusionMatrix } from '../nextflow-utils'
include { mismatchConfusionMatrix as assemblyMismatchConfusionMatrix } from '../nextflow-utils'
include { confusionMatrix as readConfusionMatrix } from '../nextflow-utils'
include { mismatchConfusionMatrix as readMismatchConfusionMatrix } from '../nextflow-utils'

// --------------------- PARAMETERS ---------------------
params.publishDir         = "results/benchmark"
params.xspectModel        = "Acinetobacter"
params.excludedSpeciesIDs = ""
params.maxForks           = 50
params.validate           = false
params.seqPlatform        = "NovaSeq"

// --------------------- WORKFLOW -----------------------
workflow {
  species_model = getModelJSON()
  name_mapping = getNameMapping(species_model)
  genomes = file("data/genomes")
  tax_report = file("data/aci_species.json")
  tax_mapping_json = strain_species_mapping(tax_report)
  assemblies = createAssemblyTable(genomes, tax_mapping_json, species_model, params.excludedSpeciesIDs)

  // Whole genome assemblies
  samples = Channel.fromPath("${genomes}/**/*.fna")
    .flatten()
  filtered_samples = assemblies
    .splitCsv(header: true, sep: '\t')
    .map { row -> row['Assembly Accession'] }
    .cross(samples.map { sample -> 
      [sample.baseName.split('_')[0..1].join('_'), sample]
    })
    .map { it[1][1] }
  classifications = classifyAssembly(filtered_samples, params.xspectModel, params.excludedSpeciesIDs)
  summarizeClassifications(assemblies, classifications.collect())
  assemblyConfusionMatrix(
    summarizeClassifications.out, 
    name_mapping, 
    'confusion_matrix.png', 
    'Xspect Acinetobacter Confusion Matrix'
  )
  assemblyMismatchConfusionMatrix(
    summarizeClassifications.out, 
    name_mapping, 
    'mismatches_confusion_matrix.png', 
    'Mismatches Confusion Matrix'
  )

  // Simulated reads
  selectForReadGen(assemblies, species_model)
  read_assemblies = selectForReadGen.out
    .splitCsv(header: true, sep: '\t')
    .map { row -> row['Assembly Accession'] }
    .cross(samples.map { sample -> 
      [sample.baseName.split('_')[0..1].join('_'), sample]
    })
    .map { it[1][1] }
  filterForChromosome(read_assemblies)
  generateReads(filterForChromosome.out)
  read_classifications = classifyRead(generateReads.out, params.xspectModel, params.excludedSpeciesIDs)
  summarizeReadClassifications(selectForReadGen.out, read_classifications.collect())
  readConfusionMatrix(
    summarizeReadClassifications.out, 
    name_mapping, 
    'read_confusion_matrix.png', 
    'Xspect Acinetobacter Read Confusion Matrix'
  )
  readMismatchConfusionMatrix(
    summarizeReadClassifications.out, 
    name_mapping, 
    'read_mismatches_confusion_matrix.png', 
    'Read Mismatches Confusion Matrix'
  )

  calculateStats(summarizeClassifications.out, summarizeReadClassifications.out)
  }

// --------------------- PROCESSES ---------------------

process getModelJSON {
  cpus 2
  memory '16 GB'

  output:
  path "species_model.json"

  script:
  """
  model_name="${params.xspectModel.toLowerCase().replaceAll('_','-')}-species.json"
  cp "\$HOME/xspect-data/models/\$model_name" species_model.json
  """
}

process getNameMapping {
  conda "conda-forge::jq"
  cpus 2
  memory '16 GB'

  input:
  path species_model

  output:
  path "name_mapping.json"

  script:
  """
  # test
  jq '.display_names | to_entries | map({key: .key, value: (.value | sub("Acinetobacter"; "A."))}) | from_entries' ${species_model} > name_mapping.json
  """

  stub:
  """
  touch name_mapping.json
  """
}


process createAssemblyTable {
  conda "conda-forge::ncbi-datasets-cli conda-forge::jq"
  cpus 2
  memory '16 GB'

  input:
  path genomes
  path tax_mapping_json
  path species_model
  val excludedSpeciesIDs

  output:
  path "assemblies.tsv"

  script:
  """
  inputfile="${genomes}/ncbi_dataset/data/assembly_data_report.jsonl"

  dataformat tsv genome --inputfile \$inputfile --fields accession,assminfo-name,organism-tax-id,assminfo-level,ani-check-status > assemblies.tsv

  # filter out assemblies with ANI check status other than "OK"
  awk -F'\t' 'NR==1 || \$5 == "OK"' assemblies.tsv > assemblies_filtered.tsv
  mv assemblies_filtered.tsv assemblies.tsv

  # add species IDs to assemblies.tsv
  declare -A species_map
  while IFS="=" read -r key val; do
    species_map["\$key"]="\$val"
  done < <(jq -r 'to_entries[] | "\\(.key)=\\(.value)"' ${tax_mapping_json})

  {
    IFS='\t' read -r -a header < assemblies.tsv
    IFS='\t'; echo -e "\${header[*]}\tSpecies ID"

    tail -n +2 assemblies.tsv | while IFS='\t' read -r acc name taxid level status; do
      species_id="\${species_map[\$taxid]:-\$taxid}"
      echo -e "\$acc\t\$name\t\$taxid\t\$level\t\$status\t\$species_id"
    done
  } > temp_assemblies.tsv
  mv temp_assemblies.tsv assemblies.tsv

  # filter out assemblies with species ID not in the species model
  jq -r '.display_names | keys | .[]' ${species_model} > valid_species.txt
  awk -F'\t' '
    BEGIN {
      while ((getline species < "valid_species.txt") > 0) {
        valid[species] = 1;
      }
      close("valid_species.txt");
    }
    NR==1 { print; next }
    \$6 in valid { print }
  ' assemblies.tsv > temp_assemblies.tsv
  mv temp_assemblies.tsv assemblies.tsv
  rm valid_species.txt

  # filter out assemblies that are part of the training set
  jq -r '.training_accessions | to_entries[] | .value[]' ${species_model} > training_accessions.txt
  awk -F'\t' '
    BEGIN {
      while ((getline acc < "training_accessions.txt") > 0) {
        training[acc] = 1;
      }
      close("training_accessions.txt");
    }
    NR==1 { print; next }
    !(\$1 in training) { print }
  ' assemblies.tsv > temp_assemblies.tsv
  mv temp_assemblies.tsv assemblies.tsv
  rm training_accessions.txt

  # filter out assemblies with excluded species IDs
  excluded_species="${excludedSpeciesIDs}"
  if [ -n "\$excluded_species" ]; then
    awk -F'\t' -v excluded="\$excluded_species" '
      BEGIN {
        # split on commas into array arr
        n = split(excluded, arr, /,/);
        for (i = 1; i <= n; i++) {
          if (arr[i] != "") {
            excluded_map[arr[i]] = 1;
          }
        }
      }
      NR==1 { print; next }
      !(\$6 in excluded_map) { print }
    ' assemblies.tsv > temp_assemblies.tsv
    mv temp_assemblies.tsv assemblies.tsv
  fi
  """

  stub:
  """
  touch assemblies.tsv
  """
}

process summarizeClassifications {
  conda "conda-forge::pandas"
  cpus 4
  memory '16 GB'
  errorStrategy 'retry'
  maxRetries 3
  publishDir params.publishDir, mode: 'copy'

  input:
  path assemblies
  path classifications

  output:
  path "classifications.tsv"

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  import json
  import os

  df = pd.read_csv('${assemblies}', sep='\\t')
  df['Prediction'] = 'unknown'

  classifications = '${classifications}'.split()

  with open(classifications[0]) as f:
    data = json.load(f)
    keys = data["scores"]["total"]
    for key in keys:
      df[str(key)] = pd.NA

  for json_file in classifications:
    basename = os.path.basename(json_file).replace('.json', '')
    accession = '_'.join(basename.split('_')[:2])
    
    with open(json_file, 'r') as f:
      data = json.load(f)
      prediction = data.get('prediction', None)

      # based on max hits if no prediction field
      if not prediction:
        hits = data.get('hits', {})
        # hits is structured as {subsequence : {species: hits}}
        total_hits = {
          species: sum(subseq.get(species, 0) for subseq in hits.values())
          for species in {s for subseq in hits.values() for s in subseq}
        }
        max_hits = max(total_hits.values())
        max_species = [species for species, species_hits in total_hits.items() if species_hits == max_hits]
        prediction = max_species[0] if len(max_species) == 1 else "ambiguous"
    
    mask = df['Assembly Accession'].str.contains(accession, na=False)
    df.loc[mask, 'Prediction'] = prediction
    
    scores = data.get('scores', {}).get('total', {})
    for species_id, score in scores.items():
      df.loc[mask, str(species_id)] = score

  df.to_csv('classifications.tsv', sep='\\t', index=False)
  """
}

process selectForReadGen {
  conda "conda-forge::pandas"
  cpus 2
  memory '16 GB'

  input:
  path assemblies
  path species_model

  output:
  path "selected_samples.tsv"

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  import json

  assemblies = pd.read_csv('${assemblies}', sep='\\t')

  training_accessions = []
  with open('${species_model}', 'r') as f:
    species_model = json.load(f)
    for id, accession in species_model["training_accessions"].items():
      training_accessions.extend(accession)
  
  assemblies = assemblies[
    (assemblies['Assembly Level'] == 'Complete Genome') |
    (assemblies['Assembly Level'] == 'Chromosome')
  ]
  assemblies = assemblies[~assemblies['Assembly Accession'].isin(training_accessions)]

  assemblies.to_csv('selected_samples.tsv', sep='\\t', index=False)
  """
}

process filterForChromosome {
  conda "bioconda::seqkit"
  cpus 2
  memory '16 GB'
  errorStrategy 'retry'
  maxRetries 3
  

  input:
  path sample

  output:
  path "${sample.baseName}_chromosome.fna"

  script:
  """
  set -euo pipefail

  seqkit sort -l -r ${sample} > sorted.tmp
  seqkit head -n 1 sorted.tmp | seqkit seq -t dna -o "${sample.baseName}_chromosome.fna"
  """
}

process generateReads {
  conda "bioconda::insilicoseq"
  cpus 4
  memory '32 GB'
  errorStrategy 'retry'
  maxRetries 3

  input:
  path sample

  output:
  path "${sample.baseName}_simulated.fq"

  script:
  """
  set -euo pipefail

  iss generate \
    --model ${params.seqPlatform} \
    --genomes "${sample}" \
    --n_reads 100000 \
    --seed 42 \
    --cpus ${task.cpus} \
    --output "${sample.baseName}_simulated"
  
  # InSilicoSeq creates paired-end files by default (_R1.fastq and _R2.fastq)
  # Concatenate them into a single file if both exist
  if [ -f "${sample.baseName}_simulated_R1.fastq" ] && [ -f "${sample.baseName}_simulated_R2.fastq" ]; then
    cat "${sample.baseName}_simulated_R1.fastq" "${sample.baseName}_simulated_R2.fastq" > "${sample.baseName}_simulated.fq"
  elif [ -f "${sample.baseName}_simulated_R1.fastq" ]; then
    mv "${sample.baseName}_simulated_R1.fastq" "${sample.baseName}_simulated.fq"
  fi
  """
}

process summarizeReadClassifications {
  conda "conda-forge::pandas"
  cpus 4
  memory '128 GB'
  errorStrategy 'retry'
  maxRetries 3
  publishDir params.publishDir, mode: 'copy'

  input:
  path read_assemblies
  path read_classifications

  output:
  path "read_classifications.tsv"

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  import json
  import os

  df_assemblies = pd.read_csv('${read_assemblies}', sep='\\t')
  
  # Create a mapping of accession to species ID
  accession_to_species = dict(zip(df_assemblies['Assembly Accession'], df_assemblies['Species ID']))
  
  classifications = '${read_classifications}'.split()
  include_header = True
  for json_file in classifications:
    basename = os.path.basename(json_file).replace('.json', '')
    accession = '_'.join(basename.split('_')[:2])
    
    species_id = accession_to_species.get(accession, 'unknown')
    
    with open(json_file, 'r') as f:
      data = json.load(f)
      scores = data.get('scores', {})
      results = []
      
      for read_name, read_scores in scores.items():
        if read_name != 'total':
          if read_scores:
            hits = data.get('hits', {}).get(read_name, {})
            max_hits = max(hits.values())
            max_species = [species for species, species_hits in hits.items() if species_hits == max_hits]
            prediction = max_species[0] if len(max_species) == 1 else "ambiguous"

            result = {
              'Assembly Accession': accession,
              'Read': read_name,
              'Prediction': prediction,
              'Species ID': species_id,
              'Rejected': True if prediction == "ambiguous" else False
            }
            
            for species, score in read_scores.items():
              result[species] = score

            results.append(result)
      
      # Reads marked as misclassified
      misclassified = data.get('misclassified', {})
      if misclassified:
        for misclass_species_id, misclass_reads in misclassified.items():
          for read_name, read_hits in misclass_reads.items():
            if read_hits:
              num_kmers = data['num_kmers'][read_name]
              read_scores = {}
              for species, hits_count in read_hits.items():
                read_scores[species] = round(hits_count / num_kmers, 2)
              
              result = {
                'Assembly Accession': accession,
                'Read': read_name,
                'Prediction': misclass_species_id,
                'Species ID': species_id,
                'Rejected': True
              }
              
              for species, score in read_scores.items():
                result[species] = score

              results.append(result)

      df_results = pd.DataFrame(results)
      df_results.to_csv('read_classifications.tsv', sep='\\t', index=False, mode='a', header=include_header)
      include_header = False
  """
}

process calculateStats {
  conda "conda-forge::pandas conda-forge::scikit-learn"
  cpus 8
  memory '256 GB'
  publishDir params.publishDir, mode: 'copy'

  input:
  path assembly_classifications
  path read_classifications

  output:
  path "stats.txt"

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  from sklearn.metrics import classification_report

  # --- Assembly ---
  df_assembly = pd.read_csv('${assembly_classifications}', sep='\\t')
  df_assembly['Species ID'] = df_assembly['Species ID'].astype(str)
  df_assembly['Prediction'] = df_assembly['Prediction'].astype(str)

  y_true_asm = df_assembly['Species ID']
  y_pred_asm = df_assembly['Prediction']

  asm_matches = (y_true_asm == y_pred_asm).sum()
  asm_total = len(df_assembly)

  asm_labels = sorted(set(y_true_asm.unique()).union(set(y_pred_asm.unique())))
  asm_report = classification_report(
      y_true_asm,
      y_pred_asm,
      labels=asm_labels,
      zero_division=0
  )

  # --- Reads ---
  df_read = pd.read_csv('${read_classifications}', sep='\\t')
  df_read['Species ID'] = df_read['Species ID'].astype(str)
  df_read['Prediction'] = df_read['Prediction'].astype(str)

  y_true_read = df_read['Species ID']
  y_pred_read = df_read['Prediction']

  read_matches = (y_true_read == y_pred_read).sum()
  read_total = len(df_read)

  read_labels = sorted(set(y_true_read.unique()).union(set(y_pred_read.unique())))
  read_report = classification_report(
      y_true_read,
      y_pred_read,
      labels=read_labels,
      zero_division=0
  )

  # --- Abstaining Metrics (Reads only) ---
  # Determine actual misclassification (prediction != ground truth)
  df_read['Actually_Misclassified'] = df_read['Species ID'] != df_read['Prediction']
  
  # Get rejection status from Rejected column
  rejected = df_read['Rejected']
  not_rejected = ~rejected
  
  # Coverage: proportion of samples that are NOT rejected
  coverage = not_rejected.sum() / read_total
  
  # Selective Accuracy: accuracy on non-rejected samples only
  if not_rejected.sum() > 0:
      selective_correct = ((df_read['Species ID'] == df_read['Prediction']) & not_rejected).sum()
      selective_accuracy = selective_correct / not_rejected.sum()
      selective_risk = 1 - selective_accuracy
  else:
      selective_accuracy = 0.0
      selective_risk = 1.0
  
  # Rejection Precision: of all rejected samples, how many were actually misclassified
  if rejected.sum() > 0:
      rejection_precision = (rejected & df_read['Actually_Misclassified']).sum() / rejected.sum()
  else:
      rejection_precision = 0.0
  
  # Rejection Recall: of all misclassified samples, how many were rejected
  if df_read['Actually_Misclassified'].sum() > 0:
      rejection_recall = (rejected & df_read['Actually_Misclassified']).sum() / df_read['Actually_Misclassified'].sum()
  else:
      rejection_recall = 0.0

  # --- Output ---
  with open('stats.txt', 'w') as f:
      f.write("=== Assembly ===\\n")
      f.write(f"Total: {asm_total}\\n")
      f.write(f"Matches: {asm_matches}\\n")
      f.write(f"Mismatches: {asm_total - asm_matches}\\n")
      f.write(f"Match Rate: {asm_matches / asm_total * 100:.2f}%\\n")
      f.write(f"Mismatch Rate: {(asm_total - asm_matches) / asm_total * 100:.2f}%\\n\\n")
      f.write("Classification report (per class):\\n")
      f.write(asm_report + "\\n")

      f.write("=== Reads ===\\n")
      f.write(f"Total: {read_total}\\n")
      f.write(f"Matches: {read_matches}\\n")
      f.write(f"Mismatches: {read_total - read_matches}\\n")
      f.write(f"Match Rate: {read_matches / read_total * 100:.2f}%\\n")
      f.write(f"Mismatch Rate: {(read_total - read_matches) / read_total * 100:.2f}%\\n\\n")
      f.write("Classification report (per class):\\n")
      f.write(read_report + "\\n")
      
      f.write("\\n=== Abstaining Metrics (Reads) ===\\n")
      f.write(f"Total Reads: {read_total}\\n")
      f.write(f"Rejected Reads: {rejected.sum()}\\n")
      f.write(f"Accepted Reads: {not_rejected.sum()}\\n")
      f.write(f"Coverage: {coverage * 100:.2f}% (proportion of non-rejected samples)\\n")
      f.write(f"Selective Accuracy: {selective_accuracy * 100:.2f}% (accuracy on non-rejected samples)\\n")
      f.write(f"Selective Risk: {selective_risk * 100:.2f}% (error rate on non-rejected samples)\\n")
      f.write(f"Rejection Precision: {rejection_precision * 100:.2f}% (of rejected, how many were truly misclassified)\\n")
      f.write(f"Rejection Recall: {rejection_recall * 100:.2f}% (of misclassified, how many were rejected)\\n")
  """
}