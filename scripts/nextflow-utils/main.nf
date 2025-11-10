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

