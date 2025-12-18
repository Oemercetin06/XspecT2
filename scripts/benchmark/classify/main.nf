process classifySample {
  conda "./scripts/nextflow-utils/environment.yml"
  cpus 4
  memory '32 GB'
  errorStrategy 'retry'
  maxRetries 3
  maxForks params.maxForks

  input:
  path sample
  val model
  val excludedSpeciesIDs

  output:
  path "${sample.baseName}.json"

  script:
  def excludeOptions = excludedSpeciesIDs ? "--exclude-species ${excludedSpeciesIDs}" : ''
  def validateFlag = params.validate ? "--validation" : ''
  """
  xspect classify species -g ${model} -i ${sample} -o ${sample.baseName}.json ${excludeOptions} ${validateFlag}
  """
}