#!/usr/bin/env nextflow

include { classifySample as classifyAssembly } from './classify'
include { classifySample as classifyRead } from './classify'

process downloadModels {
  conda "./scripts/benchmark/environment.yml"
  cpus 2
  memory '16 GB'

  output:
  path "species_model.json"

  script:
  """
  if [ ! "$HOME/xspect-data/models/acinetobacter-species.json" ]; then
    xspect models download
  fi
  cp "$HOME/xspect-data/models/acinetobacter-species.json" species_model.json
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
  path tax_report
  path species_model

  output:
  path "assemblies.tsv"

  when:
  !file("assemblies.tsv").exists()


  script:
  """
  inputfile="${genomes}/ncbi_dataset/data/assembly_data_report.jsonl"

  dataformat tsv genome --inputfile \$inputfile --fields accession,assminfo-name,organism-tax-id,assminfo-level,ani-check-status > assemblies.tsv

  # filter out assemblies with ANI check status other than "OK"
  awk -F'\t' 'NR==1 || \$5 == "OK"' assemblies.tsv > assemblies_filtered.tsv
  mv assemblies_filtered.tsv assemblies.tsv

  # map taxonmic IDs to species IDs (taxonomic IDs might be strain IDs)
  jq '
    .reports
    | map(select(.taxonomy.children != null))
    | map({
        species_id: .taxonomy.tax_id,
        children: .taxonomy.children
      })
    | map(
        . as \$entry
        | \$entry.children
        | map({ (tostring): \$entry.species_id })
        | add
      )
    | add
  ' ${tax_report} > tax_mapping.json

  # add species IDs to assemblies.tsv
  declare -A species_map
  while IFS="=" read -r key val; do
    species_map["\$key"]="\$val"
  done < <(jq -r 'to_entries[] | "\\(.key)=\\(.value)"' tax_mapping.json)

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
  publishDir "results"

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
      prediction = data.get('prediction', 'unknown')
    
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

  # use up to three assemblies for each species
  assemblies = assemblies.groupby('Species ID').head(3)

  assemblies.to_csv('selected_samples.tsv', sep='\\t', index=False)
  """
}

process generateReads {
  conda "conda-forge::pandas conda-forge::biopython"
  cpus 2
  memory '16 GB'

  input:
  path sample

  output:
  path "${sample.baseName}_simulated.fq"

  script:
  """
  #!/usr/bin/env python
  import random
  from Bio import SeqIO
  
  read_length = 100
  num_reads = 100000
  seed = 42
  
  random.seed(seed)
  sequences = list(SeqIO.parse("${sample}", "fasta"))
  chromosome_sequence = max(sequences, key=len)  # we assume the longest sequence is the chromosome
  
  ch_rec_id = chromosome_sequence.id
  ch_seq = chromosome_sequence.seq
  ch_seqlen = len(chromosome_sequence.seq)
  with open("${sample.baseName}_simulated.fq", "w") as f:
    for i in range(num_reads):
      start = random.randint(0, ch_seqlen - read_length)
      read_seq = ch_seq[start:start + read_length]
      f.write(f"@read_{i}_{ch_rec_id}_{start}-{start+read_length}\\n")
      f.write(f"{read_seq}\\n")
      f.write("+\\n")
      f.write(f"{len(read_seq)*'~'}\\n")
  """
}

process summarizeReadClassifications {
  conda "conda-forge::pandas"
  cpus 4
  memory '16 GB'
  publishDir "results"

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

  results = []
  
  classifications = '${read_classifications}'.split()
  for json_file in classifications:
    basename = os.path.basename(json_file).replace('.json', '')
    accession = '_'.join(basename.split('_')[:2])
    
    species_id = accession_to_species.get(accession, 'unknown')
    
    with open(json_file, 'r') as f:
      data = json.load(f)
      scores = data.get('scores', {})
      
      for read_name, read_scores in scores.items():
        if read_name != 'total':
          if read_scores:
            max_score = max(read_scores.values())
            max_species = [species for species, score in read_scores.items() if score == max_score]
            prediction = max_species[0] if len(max_species) == 1 else "ambiguous"

            result = {
              'Assembly Accession': accession,
              'Read': read_name,
              'Prediction': prediction,
              'Species ID': species_id
            }
            
            for species, score in read_scores.items():
              result[species] = score

            results.append(result)

  df_results = pd.DataFrame(results)
  df_results.to_csv('read_classifications.tsv', sep='\\t', index=False)
  """
}

process calculateStats {
  conda "conda-forge::pandas conda-forge::scikit-learn"
  cpus 2
  memory '16 GB'
  publishDir "results"

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
  """
}

process confusionMatrix {
  conda "conda-forge::pandas conda-forge::scikit-learn conda-forge::numpy conda-forge::matplotlib"
  cpus 2
  memory '16 GB'
  publishDir "results"

  input:
  path classifications
  path name_mapping

  output:
  path "confusion_matrix.png"

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  from sklearn.metrics import confusion_matrix
  import matplotlib.pyplot as plt
  import numpy as np
  import json
  
  df = pd.read_csv('${classifications}', sep='\\t')
  y_true = df["Species ID"].astype(str)
  y_pred = df["Prediction"].astype(str)

  with open('${name_mapping}', 'r') as f:
      name_mapping_dict = json.load(f)
  labels = list(set(y_true) | set(y_pred))
  labels = sorted(labels, key=lambda x: name_mapping_dict.get(x, x))
  display_labels = [name_mapping_dict.get(label, label) for label in labels]

  cm = confusion_matrix(y_true, y_pred, labels=labels)
  cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
  
  plt.figure(figsize=(30, 30))
  plt.imshow(cm_normalized, interpolation='nearest', cmap=plt.cm.Blues)
  plt.colorbar()
  
  plt.xticks(ticks=np.arange(len(labels)), labels=display_labels, rotation=90, fontsize=12)
  plt.yticks(ticks=np.arange(len(labels)), labels=display_labels, fontsize=12)
  
  plt.title('Xspect Acinetobacter Confusion Matrix', fontsize=24)
  plt.xlabel('Predicted Labels', fontsize=20)
  plt.ylabel('True Labels', fontsize=20)
  
  plt.savefig('confusion_matrix.png', dpi=300, bbox_inches='tight')
  """
}

process mismatchConfusionMatrix {
  conda "conda-forge::pandas conda-forge::scikit-learn conda-forge::numpy conda-forge::matplotlib"
  cpus 2
  memory '16 GB'
  publishDir "results"

  input:
  path classifications
  path name_mapping

  output:
  path "mismatches_confusion_matrix.png"

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  from sklearn.metrics import confusion_matrix
  import matplotlib.pyplot as plt
  import numpy as np
  import json

  
  df = pd.read_csv('${classifications}', sep='\\t')
  df["Species ID"] = df["Species ID"].astype(str)
  df["Prediction"] = df["Prediction"].astype(str)
  df_comparison_mismatch = df[df["Species ID"] != df["Prediction"]]

  with open('${name_mapping}', 'r') as f:
      name_mapping_dict = json.load(f)
  y_true = df_comparison_mismatch["Species ID"]
  y_pred = df_comparison_mismatch["Prediction"]
  
  labels = list(set(y_true) | set(y_pred))
  labels = sorted(labels, key=lambda x: name_mapping_dict.get(x, x))
  display_labels = [name_mapping_dict.get(label, label) for label in labels]
  
  cm = confusion_matrix(y_true, y_pred, labels=labels)
  
  plt.figure(figsize=(30, 30))
  plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
  cbar = plt.colorbar()
  cbar.ax.tick_params(labelsize=20)
  
  plt.xticks(ticks=np.arange(len(labels)), labels=display_labels, rotation=90, fontsize=16)
  plt.yticks(ticks=np.arange(len(labels)), labels=display_labels, fontsize=16)
  
  thresh = cm.max() / 2.
  for i in range(cm.shape[0]):
      for j in range(cm.shape[1]):
          plt.text(j, i, format(cm[i, j], 'd'),  # 'd' ensures integer formatting
                  horizontalalignment="center",
                  color="white" if cm[i, j] > thresh else "black",
                  fontsize=14)
  
  plt.title('Mismatches Confusion Matrix', fontsize=30)
  plt.xlabel('Predicted Labels', fontsize=24)
  plt.ylabel('True Labels', fontsize=24)
  
  plt.savefig('mismatches_confusion_matrix.png', dpi=300, bbox_inches='tight')
  """
}


workflow {
  species_model = downloadModels()
  name_mapping = getNameMapping(species_model)
  genomes = file("data/genomes")
  tax_report = file("data/aci_species.json")
  assemblies = createAssemblyTable(genomes, tax_report, species_model)

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
  classifications = classifyAssembly(filtered_samples)
  summarizeClassifications(assemblies, classifications.collect())
  confusionMatrix(summarizeClassifications.out, name_mapping)
  mismatchConfusionMatrix(summarizeClassifications.out, name_mapping)

  // Simulated reads
  selectForReadGen(assemblies, species_model)
  read_assemblies = selectForReadGen.out
    .splitCsv(header: true, sep: '\t')
    .map { row -> row['Assembly Accession'] }
    .cross(samples.map { sample -> 
      [sample.baseName.split('_')[0..1].join('_'), sample]
    })
    .map { it[1][1] }
  generateReads(read_assemblies)
  read_classifications = classifyRead(generateReads.out)
  summarizeReadClassifications(selectForReadGen.out, read_classifications.collect())

  calculateStats(summarizeClassifications.out, summarizeReadClassifications.out)
  }