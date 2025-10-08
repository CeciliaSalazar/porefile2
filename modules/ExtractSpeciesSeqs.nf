process ExtractSpeciesSeqs {

  tag "${sample_id}"
  publishDir "${params.outdir}/Species_Seqs", mode: 'copy'

  input:
  tuple val(sample_id), path(fasta_file), path(readinfo_file)

  output:
  tuple val(sample_id), path("${sample_id}.species.fa")

  script:
  """
  set -euo pipefail

  # 1) Get read IDs for species (col2 == 'S'); skip header if present
  awk -F '\\t' 'NR>1 && \$2=="S"{id=\$1; sub(/^>/,"",id); gsub(/^[ \\t]+|[ \\t]+\$/,"",id); if(id!="") print id}' \
    "!{readinfo_file}" > species.ids

  # If none, create empty fasta and exit cleanly
  if ! grep -qve '^[[:space:]]*\$' species.ids 2>/dev/null; then
    : > "!{sample_id}.species.fa"
    echo "[ExtractSpeciesSeqs] sample=!{sample_id} species_ids=0 sequences_written=0" >&2
    exit 0
  fi

  # 2) Filter FASTA (normalize both sides)
  awk '
    function norm(s,t){
      t=s
      gsub(/^[ \\t]+|[ \\t]+\$/,"",t)
      sub(/^>/,"",t)
      sub(/[ \\t].*\\$/,"",t)
      sub(/\\|.*\\$/,"",t)
      sub(/\\/.*\\$/,"",t)
      sub(/\\.[0-9]+\\$/,"",t)
      return t
    }
    BEGIN{
      while((getline line < "species.ids")>0){
        if(line ~ /^[ \\t]*\$/) continue
        id = norm(line)
        if(id!="") ids[id]=1
      }
      close("species.ids")
    }
    /^>/{
      hdr = substr(\$0,2)
      keep = (norm(hdr) in ids)
    }
    { if(keep) print \$0 }
  ' "!{fasta_file}" > "!{sample_id}.species.fa"

  n_ids=\$(grep -cvE '^[[:space:]]*\$' species.ids || true)
  n_seq=\$(grep -c '^>' "!{sample_id}.species.fa" || true)
  echo "[ExtractSpeciesSeqs] sample=!{sample_id} species_ids=\${n_ids} sequences_written=\${n_seq}" >&2
  """
}
