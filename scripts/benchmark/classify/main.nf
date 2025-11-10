process classifySample {
  conda "./scripts/benchmark/environment.yml"
  cpus 4
  memory '32 GB'

  input:
  path sample
  val model

  output:
  path "${sample.baseName}.json"

  script:
  """
  xspect classify species -g ${model} -i ${sample} -o ${sample.baseName}.json
  """

  stub:
  """
  mkdir -p results
  touch results/${sample.baseName}.json
  """
}