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

  # 1) Leer IDs con flag [S] del .read_info (tab-delimited; 1ª columna = read id)
  #    Evitamos clases POSIX y usamos patrones simples para recortar espacios.
  awk -F '\t' 'NR==1{next} /\[S\]/{
      id=$1
      sub(/^>/,"",id)
      sub(/^[ \t\r\n]+/,"",id)
      sub(/[ \t\r\n]+$/,"",id)
      if(id!="") print id
  }' "!{readinfo_file}" > species.ids

  # Si no hay IDs, crear salida vacía y salir limpio
  if ! grep -qve '^[ \t\r\n]*$' species.ids 2>/dev/null; then
    : > "!{sample_id}.species.fa"
    exit 0
  fi

  # 2) Filtrar FASTA por IDs normalizando cabeceras (antes de espacio, '|' o '/')
  awk '
    function norm(s, t){
      t=s
      sub(/[ \t].*$/, "", t)
      sub(/\|.*$/,   "", t)
      sub(/\/.*$/,   "", t)
      return t
    }
    BEGIN{
      while((getline line < "species.ids")>0){
        if(line ~ /^[ \t\r\n]*$/) continue
        id = norm(line)
        ids[id]=1
      }
      close("species.ids")
      keep=0
    }
    /^>/{
      hdr  = substr($0,2)
      keep = (norm(hdr) in ids)
    }
    { if(keep) print $0 }
  ' "!{fasta_file}" > "!{sample_id}.species.fa"

  # Log
  n_ids=$(grep -cvE '^[ \t\r\n]*$' species.ids || true)
  n_seq=$(grep -c '^>' "!{sample_id}.species.fa" || true)
  echo "[ExtractSpeciesSeqs] sample=!{sample_id} species_ids=$n_ids sequences_written=$n_seq" >&2
  '''
}
