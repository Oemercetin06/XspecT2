#!/usr/bin/env nextflow

process strain_species_mapping {
    conda "conda-forge::jq"

    input:
    path tax_report

    output:
    path "tax_mapping.json", emit: 'tax_mapping_json'

    script:
    """
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
  """
}

process confusionMatrix {
  conda "conda-forge::pandas conda-forge::scikit-learn conda-forge::numpy conda-forge::matplotlib"
  cpus 4
  memory '32 GB'
  publishDir "${params.publishDir ?: 'results'}", mode: 'copy', enabled: params.publishDir != null

  input:
  path classifications
  path name_mapping
  val output_filename
  val title

  output:
  path "${output_filename}"

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  from sklearn.metrics import confusion_matrix
  import matplotlib.pyplot as plt
  import numpy as np
  import json
  
  df = pd.read_csv('${classifications}', sep='\\t', usecols=['Species ID', 'Prediction'], dtype=str)
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
  
  plt.title('${title}', fontsize=24)
  plt.xlabel('Predicted Species', fontsize=20)
  plt.ylabel('NCBI-annotated Species', fontsize=20)
  
  plt.savefig('${output_filename}', dpi=300, bbox_inches='tight')
  """
}

process mismatchConfusionMatrix {
  conda "conda-forge::pandas conda-forge::scikit-learn conda-forge::numpy conda-forge::matplotlib"
  cpus 4
  memory '32 GB'
  publishDir "${params.publishDir ?: 'results'}", mode: 'copy', enabled: params.publishDir != null

  input:
  path classifications
  path name_mapping
  val output_filename
  val title

  output:
  path "${output_filename}"

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  from sklearn.metrics import confusion_matrix
  import matplotlib.pyplot as plt
  import numpy as np
  import json

  
  df = pd.read_csv('${classifications}', sep='\\t', usecols=['Species ID', 'Prediction'], dtype=str)

  df_comparison_mismatch = df[df["Species ID"] != df["Prediction"]]
  if df_comparison_mismatch.empty:
      print("No mismatches found. Skipping mismatch confusion matrix generation.")
      with open('${output_filename}', 'w') as f:
          f.write('')
      exit(0)

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
  
  plt.title('${title}', fontsize=30)
  plt.xlabel('Predicted Species', fontsize=24)
  plt.ylabel('NCBI-annotated Species', fontsize=24)
  
  plt.savefig('${output_filename}', dpi=300, bbox_inches='tight')
  """
}

