process classifySample {
  conda "./scripts/benchmark/environment.yml"
  cpus 4
  memory '32 GB'

  input:
  path sample

  output:
  path "${sample.baseName}.json"

  script:
  """
  xspect classify species -g Acinetobacter -i ${sample} -o ${sample.baseName}.json
  """

  stub:
  """
  mkdir -p results
  touch results/${sample.baseName}.json
  """
}