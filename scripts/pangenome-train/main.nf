#!/usr/bin/env nextflow

include { strain_species_mapping } from '../nextflow-utils'

// --------------------- PARAMETERS ---------------------
params.dataset_dir = "data/genomes/ncbi_dataset"
params.cpus = 32
// Clustering method: either "98" or "GF"
params.clusteringMethod = "98"

// --------------------- WORKFLOW -----------------------
workflow {
  // Input validation
  if (!['98', 'GF'].contains(params.clusteringMethod)) {
    throw new IllegalArgumentException("Invalid params.clusteringMethod: '${params.clusteringMethod}'. Must be '98' for 98% similarity clustering or 'GF' for gene family clustering.")
  }

  // Build maps
  taxon_map = EXTRACT_TAXON(file(params.dataset_dir) + "/data/assembly_data_report.jsonl")
  fasta_index = INDEX_FASTAS(file(params.dataset_dir))
  map = JOIN_MAPS(taxon_map, fasta_index)

  tax_report = file("data/aci_species.json")

  tax_mapping_json = strain_species_mapping(tax_report)
  reassigned_map = REASSIGN_MAP(map, tax_mapping_json)

  // Exclude taxa with ambiguous or provisional names (e.g., containing "sp." or "Candidatus")
  filtered_name_map = FILTER_TAXON_NAMES(reassigned_map, tax_report)
  // Filter to only ACB clade species
  // acb_map = FILTER_ACB(reassigned_map)

  split_files = TRAIN_TEST_SPLIT(filtered_name_map)
  filtered = FILTER_SINGLETON_TAXA(split_files.train.flatten().collect())
  train_taxon_files = filtered.train.flatten()

  pangenomes = PANGENOME(train_taxon_files, file(params.dataset_dir))
  pangenome_fastas = EXTRACT_FASTA(pangenomes)
  if (params.clusteringMethod == '98') {
    clustered_fastas = CLUSTER_FASTAS(pangenome_fastas)
  } else {
    clustered_fastas = pangenome_fastas
  }
  organized_fasta_folders = ORGANIZE_PANGENOME_FASTAS(clustered_fastas.collect())

  fasta_folder_ch = organized_fasta_folders.core.mix(
    organized_fasta_folders.softcore_95,
    organized_fasta_folders.softcore_90,
    organized_fasta_folders.softcore_85,
    organized_fasta_folders.softcore_80,
    organized_fasta_folders.softcore_75,
    organized_fasta_folders.persistent,
    organized_fasta_folders.persistent_shell,
    organized_fasta_folders.persistent_shell_cloud,
  )

  trained_models = XSPECT_TRAIN(fasta_folder_ch)
  UPDATE_MODEL_METADATA(trained_models.collect(), train_taxon_files.collect(), tax_report)
}

// --------------------- PROCESSES ----------------------

process EXTRACT_TAXON {
  publishDir "results/pangenome-train/maps", mode: 'copy'

  conda "conda-forge::jq"

  input:
  path jsonl

  output:
  path "genome_to_taxid.tsv", emit: 'genome_to_taxid_tsv'

  script:
  """
    set -euo pipefail
    jq -r 'select(.averageNucleotideIdentity.taxonomyCheckStatus == "OK" and .checkmInfo.completeness >= 95) | [.accession, (.organism.taxId // .taxon.taxId // "Unknown")] | @tsv' "${jsonl}" > genome_to_taxid.tsv
    """
}

process INDEX_FASTAS {
  publishDir "results/pangenome-train/maps", mode: 'copy'

  input:
  path dataset_dir

  output:
  path "accession_to_fna.tsv", emit: 'accession_to_fna_tsv'

  script:
  """
    set -euo pipefail
    find "${dataset_dir}/data" -type f '(' -name "*genomic.fna" -o -name "*genomic.fna.gz" ')' |
      awk '{
        n = split(\$0, a, "/");
        acc = a[n-1];
        print acc "\t" \$0
      }' |
      sort > accession_to_fna.tsv
    """
}

process JOIN_MAPS {
  publishDir "results/pangenome-train/maps", mode: 'copy'

  input:
  path genome_to_taxid_tsv
  path accession_to_fna_tsv

  output:
  path "accession_taxid_fna.tsv", emit: 'accession_taxid_fna_tsv'

  script:
  """
    set -euo pipefail
    join -t \$'\\t' -1 1 -2 1 <(sort "${genome_to_taxid_tsv}") <(sort "${accession_to_fna_tsv}") > accession_taxid_fna.tsv
    """
}

process REASSIGN_MAP {
  publishDir "results/pangenome-train/maps", mode: 'copy'
  conda "conda-forge::jq"

  input:
  path accession_taxid_fna_tsv
  path tax_mapping_json

  output:
  path "accession_taxid_fna_reassigned.tsv", emit: 'accession_taxid_fna_reassigned_tsv'

  script:
  """
    set -euo pipefail
    awk -F'\\t' 'NR==FNR{map[\$1]=\$2; next}{
        if (\$2 in map) { \$2 = map[\$2] }
        print \$0
      }' <(jq -r 'to_entries[] | "\\(.key)\\t\\(.value)"' "${tax_mapping_json}") "${accession_taxid_fna_tsv}" > accession_taxid_fna_reassigned.tsv
    """
}

process FILTER_TAXON_NAMES {
  publishDir "results/pangenome-train/maps", mode: 'copy'
  conda "conda-forge::jq"

  input:
  path reassigned_tsv
  path tax_report

  output:
  path "accession_taxid_fna_filtered.tsv", emit: 'accession_taxid_fna_filtered_tsv'

  script:
  """
    set -euo pipefail

    # Extract valid species-level tax IDs whose scientific names do NOT contain 'sp.' or 'Candidatus' (case-insensitive)
    jq -r '.reports[]
      | .taxonomy
      | select((.current_scientific_name.name // "") | test("(?i)\\\\bCandidatus\\\\b|\\\\bsp\\\\.", "i") | not)
      | .tax_id' "${tax_report}" \
      | sort -u > valid_species_taxids.txt

    # Keep only rows where the (possibly reassigned) taxid is in the valid list
    awk -F'\t' 'NR==FNR{valid[\$1]=1; next} (\$2 in valid)' valid_species_taxids.txt "${reassigned_tsv}" \
      > accession_taxid_fna_filtered.tsv
  """
}

process FILTER_ACB {
  publishDir "results/pangenome-train/maps", mode: 'copy'

  input:
  path reassigned_tsv

  output:
  path "accession_taxid_fna_ACB.tsv", emit: 'accession_taxid_fna_acb_tsv'

  script:
  """
    set -euo pipefail
    awk -F'\t' '(\$2=="470"||\$2=="48296"||\$2=="106654"||\$2=="471"||\$2=="1530123"||\$2=="1785128")' "${reassigned_tsv}" > accession_taxid_fna_ACB.tsv
    """
}

process TRAIN_TEST_SPLIT {
  publishDir "results/pangenome-train/maps/by_taxid", mode: 'copy'

  input:
  path accession_taxid_fna_tsv

  output:
  path "*_train.tsv", emit: 'train'
  path "*_test.tsv", emit: 'test'

  script:
  """
    #!/usr/bin/env python
    import random
    from collections import defaultdict

    # Group rows by taxid
    groups = defaultdict(list)
    with open("${accession_taxid_fna_tsv}", 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\\t', 2)
            if len(parts) < 3:
                continue
            acc = parts[0]
            tax = parts[1]
            path = parts[2]
            groups[tax].append(f"{acc}\\t{path}")

    # For each taxid, shuffle and split into train (80%) and test (20%)
    random.seed(42)
    for tax, rows in groups.items():
        random.shuffle(rows)
        
        split_index = max(1, (8 * len(rows)) // 10)
        train_rows = rows[:split_index]
        test_rows = rows[split_index:]
        
        with open(f"{tax}_train.tsv", 'w') as f:
            f.write('\\n'.join(train_rows) + '\\n')
        
        with open(f"{tax}_test.tsv", 'w') as f:
            f.write('\\n'.join(test_rows) + '\\n')
    """
}


process FILTER_SINGLETON_TAXA {
  publishDir "results/pangenome-train/maps/by_taxid", mode: 'copy', pattern: "*_train.tsv"
  publishDir "results/pangenome-train/maps", mode: 'copy', pattern: "singletons/*.tsv"

  input:
  path train_file

  output:
  path "filtered_train/*_train.tsv", emit: 'train'
  path "singletons/*.tsv", emit: 'singletons'

  script:
  """
  #!/usr/bin/env python3
  import os
  import shutil

  train_files = "${train_file}".split()
  for train_file in train_files:
    base = os.path.basename(train_file)
    taxid = base.split('_', 1)[0]

    # Count non-empty lines
    n = 0
    with open(train_file, 'r') as f:
      for line in f:
        if line.strip():
          n += 1

    if n > 3:
      os.makedirs('filtered_train', exist_ok=True)
      shutil.copy(train_file, os.path.join('filtered_train', base))
    else:
      os.makedirs('singletons', exist_ok=True)
      shutil.copy(train_file, os.path.join('singletons', f"{taxid}.tsv"))
  """
}


process PANGENOME {
  publishDir "results/pangenome-train/pangenome", mode: 'copy'
  conda "bioconda::ppanggolin"
  cpus { genomes_list.baseName == "470_train" ? params.cpus * 2 : params.cpus }
  memory params.cpus * 4 + " GB"

  input:
  path genomes_list
  path dataset_dir

  output:
  path "${genomes_list.baseName.split('_')[0]}.h5", emit: 'pangenome_h5'

  script:
  """
    set -euo pipefail
    ppanggolin all --fasta "${genomes_list}" --cpu ${task.cpus}
    mv **/*.h5 ${genomes_list.baseName.split('_')[0]}.h5
    """
}

process EXTRACT_FASTA {
  publishDir "results/pangenome-train/${params.clusteringMethod}-fastas", mode: 'copy'
  conda "bioconda::ppanggolin"
  cpus params.cpus
  memory params.cpus * 4 + " GB"

  input:
  path pangenome_h5

  output:
  path "${pangenome_h5.baseName}"


  script:
  def geneOpt = params.clusteringMethod == '98' ? '--genes' : '--gene_families'
  """
    set -euo pipefail

    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}" ${geneOpt} core -f
    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}" ${geneOpt} softcore -f
    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}/sc90" ${geneOpt} softcore -f --soft_core 0.9
    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}/sc85" ${geneOpt} softcore -f --soft_core 0.85
    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}/sc80" ${geneOpt} softcore -f --soft_core 0.8
    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}/sc75" ${geneOpt} softcore -f --soft_core 0.75
    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}" ${geneOpt} persistent -f
    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}" ${geneOpt} shell -f
    ppanggolin fasta -p "${pangenome_h5}" --output "${pangenome_h5.baseName}" ${geneOpt} cloud -f
    """
}

process CLUSTER_FASTAS {
  publishDir "results/pangenome-train/clustered", mode: 'copy'
  conda "bioconda::cd-hit bioconda::seqkit"
  cpus params.cpus
  memory params.cpus * 4 + " GB"

  input:
  path fasta_dir

  output:
  path "${fasta_dir.baseName}_clustered"

  script:
  """
    set -euo pipefail
    
    mkdir -p "${fasta_dir.baseName}_clustered"
    
    find -L "${fasta_dir}" -type f -name "*.fna" | while read -r fasta_file; do
      rel_path=\${fasta_file#${fasta_dir}/}
      output_file="${fasta_dir.baseName}_clustered/\${rel_path}"
      temp_file="\${output_file}.temp"
      
      mkdir -p "\$(dirname "\${output_file}")"
      
      seqkit rmdup -s "\${fasta_file}" -o "\${temp_file}"
      
      # Cluster using cd-hit-est
      # -c 0.98: 98% sequence identity threshold
      # -T 0: use all available threads
      # -M 0: unlimited memory
      # -d 0: keep full sequence descriptions
      cd-hit-est -i "\${temp_file}" -o "\${output_file}" \\
        -c 0.98 -T ${task.cpus} -M 0 -d 0
      
      rm -f "\${temp_file}" "\${output_file}.clstr"
    done
    """
}

process ORGANIZE_PANGENOME_FASTAS {
  publishDir "results/pangenome-train/${params.clusteringMethod}-organized", mode: 'copy'

  input:
  path pangenome_fastas

  output:
  path "core", emit: core
  path "softcore_95", emit: softcore_95
  path "softcore_90", emit: softcore_90
  path "softcore_85", emit: softcore_85
  path "softcore_80", emit: softcore_80
  path "softcore_75", emit: softcore_75
  path "persistent", emit: persistent
  path "persistent_shell", emit: persistent_shell
  path "persistent_shell_cloud", emit: persistent_shell_cloud

  script:
  """
    #!/usr/bin/env python
    import os
    import shutil
    import glob
    import re

    pangenome_dirs = "${pangenome_fastas}".split()
    
    for fasta_dir in pangenome_dirs:
        if not os.path.isdir(fasta_dir):
            continue
            
        taxid = os.path.basename(fasta_dir).split('_')[0]
        
        # Find all fasta files in the directory and subdirectories
        files = []
        files.extend(glob.glob(f"{fasta_dir}/**/*.fna", recursive=True))
        files.extend(glob.glob(f"{fasta_dir}/**/*.fasta", recursive=True))
        for fasta in sorted(set(files)):
            base = os.path.basename(fasta)
            
            # Determine partition type from filename
            out_dirs = []
            if "core" in base and "softcore" not in base:
                out_dirs.append("core")
            elif "softcore" in base:
                # Determine softcore threshold from path or filename
                # Default softcore (0.95) - check if it's in the main directory
                parent_dir = os.path.basename(os.path.dirname(fasta))
                
                if parent_dir == "sc90":
                    out_dirs.append("softcore_90")
                elif parent_dir == "sc85":
                    out_dirs.append("softcore_85")
                elif parent_dir == "sc80":
                    out_dirs.append("softcore_80")
                elif parent_dir == "sc75":
                    out_dirs.append("softcore_75")
                else:
                    # Default softcore (0.95)
                    out_dirs.append("softcore_95")
            elif "persistent" in base:
                out_dirs.append("persistent")
                out_dirs.append("persistent_shell")
                out_dirs.append("persistent_shell_cloud")
            elif "shell" in base:
                out_dirs.append("persistent_shell")
                out_dirs.append("persistent_shell_cloud")
            elif "cloud" in base:
                out_dirs.append("persistent_shell_cloud")

            for out_dir in out_dirs:
                dest_dir = f"{out_dir}/cobs/{taxid}"
                os.makedirs(dest_dir, exist_ok=True)
                dest_path = os.path.join(dest_dir, base)
                shutil.copy(fasta, dest_path)
    """
}

process XSPECT_TRAIN {
  conda "./scripts/nextflow-utils/environment.yml"
  cpus params.cpus
  memory params.cpus * 4 + " GB"

  input:
  path organized_fasta_folder

  output:
  path "model_${organized_fasta_folder.baseName}.txt"

  script:
  """
    set -euo pipefail

    xspect models train directory -g "Acinetobacter_${organized_fasta_folder.baseName}_${params.clusteringMethod}" -i "${organized_fasta_folder}"
    echo "Acinetobacter_${organized_fasta_folder.baseName}_${params.clusteringMethod}" > "model_${organized_fasta_folder.baseName}.txt"
    """
}

process UPDATE_MODEL_METADATA {
  publishDir "results/pangenome-train/metadata", mode: 'copy'

  input:
  path model_names
  path train_files
  path tax_report

  output:
  path "training_accessions.json"

  script:
  """
    #!/usr/bin/env python
    import json
    import os
    from pathlib import Path

    # Load taxonomy data to get scientific names
    with open("${tax_report}", 'r') as f:
        tax_data = json.load(f)
    
    # Build a mapping from tax_id to scientific name
    taxid_to_name = {}
    for report in tax_data.get('reports', []):
        taxonomy = report.get('taxonomy', {})
        tax_id = taxonomy.get('tax_id')
        scientific_name = taxonomy.get('current_scientific_name', {}).get('name')
        if tax_id and scientific_name:
            taxid_to_name[str(tax_id)] = scientific_name

    # Parse all training files and collect accessions by taxid
    training_data = {}
    train_files = "${train_files}".split()
    
    for train_file in train_files:
        if not os.path.exists(train_file):
            continue
        
        # Extract taxid from filename (e.g., "470_train.tsv")
        taxid = Path(train_file).stem.split('_')[0]
        
        # Read accessions from the training file
        accessions = []
        with open(train_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\\t')
                    if parts:
                        accessions.append(parts[0])
        
        training_data[taxid] = accessions
    
    # Save the training accessions mapping
    with open('training_accessions.json', 'w') as f:
        json.dump(training_data, f, indent=4)
    
    # Update model metadata files
    model_names = "${model_names}".split()
    for model_name_file in model_names:
        with open(model_name_file, 'r') as f:
            model_name = f.read().strip()
        
        # Determine model path based on naming convention
        # Extract the model suffix from model name (e.g., "Acinetobacter_core_98" or "Acinetobacter_core_GF")
        model_suffix = model_name.split('_', maxsplit=1)[-1].replace('_', '-').lower()
        
        # Construct model file path
        model_file = Path.home() / "xspect-data" / "models" / f"acinetobacter-{model_suffix}-species.json"
        
        if model_file.exists():
            # Read existing model metadata
            with open(model_file, 'r') as f:
                metadata = json.load(f)
            
            # Update with training accessions
            metadata['training_accessions'] = training_data
            
            # Update display_names with tax_id: scientific_name mappings
            display_names = {}
            for taxid in training_data.keys():
                if taxid in taxid_to_name:
                    display_names[taxid] = taxid_to_name[taxid]
            metadata['display_names'] = display_names
            
            # Write back to model file
            with open(model_file, 'w') as f:
                json.dump(metadata, f, indent=4)
            
            print(f"Updated {model_file}")
        else:
            print(f"Model file not found: {model_file}")
    """
}
