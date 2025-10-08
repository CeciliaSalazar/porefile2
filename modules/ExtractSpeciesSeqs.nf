process ExtractSpeciesSeqs {

  input:
    tuple val(sample_id), path(fasta_file), path(readinfo_file)

  output:
    tuple val(sample_id), path("${sample_id}.species.fa")

  script:
  """
  set -euo pipefail

  # 1) IDs a nivel especie: .read_info es TAB-delimited; col 1 = read_id; alguna columna contiene '[S]'
  awk -F '\\t' 'NR==1{next} /\\[S\\]/{id=\\$1; sub(/^>/,"",id); gsub(/^[[:space:]]+/, "", id); gsub(/[[:space:]]+$/, "", id); if(id!="") print id}' "!{readinfo_file}" > species.ids

  # Si no hay IDs, crear salida vacía y salir limpio
  if ! grep -qve '^[[:space:]]*$' species.ids 2>/dev/null; then
    : > "!{sample_id}.species.fa"
    exit 0
  fi

  # 2) Filtrar FASTA por esos IDs normalizando el header (corta en 1er espacio, '|' o '/')
  awk '
    function norm(s,t){
      t=s
      sub(/[ \\t].*$/, "", t)   # corta en primer espacio
      sub(/\\|.*$/,   "", t)    # corta en primer |
      sub(/\\/.*$/,   "", t)    # corta en primer /
      return t
    }
    BEGIN{
      while((getline line < "species.ids")>0){
        if(line ~ /^[[:space:]]*$/) continue
        id = norm(line)
        ids[id]=1
      }
      close("species.ids")
      keep=0
    }
    /^>/{
      hdr  = substr(\\$0,2)
      keep = (norm(hdr) in ids)
    }
    { if(keep) print \\$0 }
  ' "!{fasta_file}" > "!{sample_id}.species.fa"

  # Log
  n_ids=\$(grep -cvE '^[[:space:]]*$' species.ids || true)
  n_seq=\$(grep -c '^>' "!{sample_id}.species.fa" || true)
  echo "[ExtractSpeciesSeqs] sample=!{sample_id} species_ids=\$n_ids sequences_written=\$n_seq" >&2
  """
}

