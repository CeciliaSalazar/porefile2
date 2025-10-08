process ExtractSpeciesSeqs {

  tag "${sample_id}"
  publishDir "${params.outdir}/Species_Seqs", mode: 'copy'

  input:
  tuple val(sample_id), path(fasta_file), path(readinfo_file)

  output:
  tuple val(sample_id), path("${sample_id}.species.fa")

  shell:
  '''
  set -euo pipefail
  # Minimal no-op to prove parsing works
  : > "!{sample_id}.species.fa"
  '''
}
