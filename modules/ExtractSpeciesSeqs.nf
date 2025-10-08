process ExtractSpeciesSeqs {

  tag "${sample_id}"
  publishDir "${params.outdir}/Species_Seqs", mode: 'copy'

  input:
  tuple val(sample_id), path(fasta_file), path(readinfo_file)

  output:
  tuple val(sample_id), path("${sample_id}.species.fa")

  /*
   * Extract read IDs at species level from .read_info and subset the FASTA.
   * .read_info (TAB-delimited):
   *   col1 = read_id
   *   col2 = rank letter (S/G/F/...)
   *   later columns contain taxonomy and may include "[S] Genus species"
   */
  script:
  """
  set -euo pipefail

  # 1) Collect species-level IDs: rank letter 'S' OR line contains "[S]"
  awk -f /dev/stdin "${readinfo_file}" > species.ids <<'AWK'
  BEGIN { FS = "\t"; OFS = "\t" }
  NR==1 { next }
  ($2 == "S") || (index($0,"[S]") > 0) {
    id = $1
    sub(/^>/, "", id)
    gsub(/^[ \t]+|[ \t]+$/, "", id)
    if (id != "") print id
  }
AWK

  # If no IDs, emit empty fasta and exit cleanly
  if ! grep -qve '^[[:space:]]*$' species.ids 2>/dev/null; then
    : > "${sample_id}.species.fa"
    echo "[ExtractSpeciesSeqs] sample=${sample_id} species_ids=0 sequences_written=0 (no species-level)" >&2
    exit 0
  fi

  # 2) Filter FASTA by normalized headers
  awk -f /dev/stdin "${fasta_file}" > "${sample_id}.species.fa" <<'AWK'
  function norm_once(s) {
    gsub(/^[ \t]+|[ \t]+$/, "", s)
    sub(/^>/, "", s)
    sub(/[ \t].*$/, "", s)
    sub(/\|.*$/, "", s)
    sub(/\/.*$/, "", s)
    sub(/\.[0-9]+$/, "", s)
    return s
  }
  function norm(s) { return norm_once(norm_once(s)) }

  BEGIN{
    while ( (getline line < "species.ids") > 0 ) {
      if (line ~ /^[ \t]*$/) continue
      id = norm(line)
      if (id != "") ids[id] = 1
    }
    close("species.ids")
    keep = 0
  }
  /^>/{
    hdr = substr($0,2)
    keep = (norm(hdr) in ids)
  }
  { if (keep) print $0 }
AWK

  n_ids=$(grep -cvE '^[[:space:]]*$' species.ids || true)
  n_seq=$(grep -c '^>' "${sample_id}.species.fa" || true)
  echo "[ExtractSpeciesSeqs] sample=${sample_id} species_ids=${n_ids} sequences_written=${n_seq}" >&2
  """
}
