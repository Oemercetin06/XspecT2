#!/usr/bin/env nextflow

// --------------------- PARAMETERS ---------------------
params.classifications = "results/benchmark/read_classifications.tsv"
params.publishDir      = "results/score-svm"
params.test_size       = 0.2
params.random_state    = 42
params.sample_size     = 100000
params.mismatch_fraction = 0.2

// --------------------- WORKFLOW -----------------------
workflow {
  classifications_file = file(params.classifications)

  DOWNSAMPLE(classifications_file)
  UNSORTED_SVM(DOWNSAMPLE.out.downsampled)
}

// --------------------- PROCESSES ----------------------

process DOWNSAMPLE {
  conda "conda-forge::pandas"
  cpus 16
  memory '256 GB'

  input:
    path classifications

  output:
    path "downsampled.tsv", emit: downsampled

  script:
  """
  #!/usr/bin/env python
  import pandas as pd

  df = pd.read_csv("${classifications}", sep="\\t")

  # apply to_numeric to reduce memory usage
  score_cols = [col for col in df.columns if col.isdigit()]
  for col in score_cols:
      df[col] = pd.to_numeric(df[col], downcast='float')

  # Same rows may already be rejected due to ambiguous classifications - filter them out
  df = df[df["Rejected"] == False].copy()
      
  df["label"] = (df["Species ID"].astype(str) == df["Prediction"].astype(str)).astype(int)
  
  sample_size = ${params.sample_size}
  mismatch_fraction = ${params.mismatch_fraction}
  random_state = ${params.random_state}
  
  # Separate mismatches (label=0) and correct predictions (label=1)
  mismatches = df[df["label"] == 0].copy()
  correct = df[df["label"] == 1].copy()
  
  # Calculate target sizes for each class
  n_mismatches = int(sample_size * mismatch_fraction)
  n_correct = sample_size - n_mismatches
  
  print(f"Target: {n_mismatches} mismatches, {n_correct} correct predictions")
  print(f"Available: {len(mismatches)} mismatches, {len(correct)} correct predictions")
  
  # Function to sample with even species distribution using groupby
  def stratified_sample_by_species(data, n_samples, random_state):
      if len(data) <= n_samples:
          return data
      
      n_species = data["Species ID"].nunique()
      samples_per_species = n_samples // n_species
      
      # Sample evenly from each species
      sampled = data.groupby("Species ID", group_keys=False).apply(
          lambda x: x.sample(n=min(samples_per_species, len(x)), random_state=random_state)
      )
      
      # If we need more samples, sample remaining from the full dataset
      if len(sampled) < n_samples:
          remaining_needed = n_samples - len(sampled)
          additional = data.drop(sampled.index).sample(n=min(remaining_needed, len(data) - len(sampled)), random_state=random_state)
          sampled = pd.concat([sampled, additional])
      
      return sampled.sample(frac=1, random_state=random_state).reset_index(drop=True)
  
  # Sample each class with stratification by species
  sampled_mismatches = stratified_sample_by_species(mismatches, n_mismatches, random_state)
  sampled_correct = stratified_sample_by_species(correct, n_correct, random_state + 1000)
  
  # Combine and shuffle
  df_sampled = pd.concat([sampled_mismatches, sampled_correct], ignore_index=True)
  df_sampled = df_sampled.sample(frac=1, random_state=random_state).reset_index(drop=True)
  
  print(f"Final dataset: {len(df_sampled)} samples")
  print(f"  Mismatches: {(df_sampled['label'] == 0).sum()} ({(df_sampled['label'] == 0).sum() / len(df_sampled):.2%})")
  print(f"  Correct: {(df_sampled['label'] == 1).sum()} ({(df_sampled['label'] == 1).sum() / len(df_sampled):.2%})")
  print(f"  Unique species in mismatches: {df_sampled[df_sampled['label'] == 0]['Species ID'].nunique()}")
  print(f"  Unique species in correct: {df_sampled[df_sampled['label'] == 1]['Species ID'].nunique()}")
  
  df_sampled.to_csv("downsampled.tsv", sep="\\t", index=False)
  """
}

process UNSORTED_SVM {
  publishDir params.publishDir, mode: 'copy'
  conda "conda-forge::python=3.12 conda-forge::pandas conda-forge::scikit-learn conda-forge::numpy"
  cpus 16
  memory '256 GB'

  input:
    path classifications

  output:
    path "classification_report.txt", emit: classification_report

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  from sklearn.svm import SVC
  from sklearn.model_selection import train_test_split
  from sklearn.metrics import classification_report

  df = pd.read_csv("${classifications}", sep="\\t")

  # apply use_to_numeric to reduce memory usage
  score_cols = [col for col in df.columns if col.isdigit()]
  for col in score_cols:
      df[col] = pd.to_numeric(df[col], downcast='float')
      
  df["label"] = (df["Species ID"].astype(str) == df["Prediction"].astype(str)).astype(int)
  
  feature_cols = [c for c in df.columns if c.isdigit()]
  X = df[feature_cols].values
  y = df["label"].values
  
  X_train, X_test, y_train, y_test = train_test_split(
      X, y, test_size=${params.test_size}, random_state=${params.random_state}, stratify=y
  )

  svm = SVC(kernel="rbf", probability=True, class_weight="balanced", random_state=42)
  svm.fit(X_train, y_train)
  
  y_pred = svm.predict(X_test)
  
  with open('classification_report.txt', 'w') as f:
    f.write(f"SVM Classification Report\\n")
    f.write(f"========================\\n\\n")
    f.write(f"Total samples used: {len(df)}\\n")
    f.write(f"Training samples: {len(X_train)}\\n")
    f.write(f"Test samples: {len(X_test)}\\n\\n")
    f.write(classification_report(y_test, y_pred, target_names=["Incorrect", "Correct"]))
  """
}

