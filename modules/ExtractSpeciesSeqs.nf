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

  # 1) Extraer IDs a nivel especie desde .read_info (TAB-delimited)
  #    - Col 1 = read_id
  #    - Col 2 = rank (letra); especie si $2 == "S"
  #    - O si la línea contiene "[S]" en cualquier parte
  #    - Saltar cabecera (NR==1)
  awk -F '\t' '
    NR==1 { next }
    ($2=="S") || (index($0,"[S]")>0) {
      id=$1
      sub(/^>/,"",id)
      gsub(/^[ \t]+|[ \t]+$/,"",id)
      if (id!="") print id
    }
  ' "!{readinfo_file}" > species.ids

  # Copia para depurar (opcional)
  cp species.ids "!{sample_id}.species.ids" || true

  # Si no hay IDs, salida vacía y terminar limpio
  if ! grep -qve '^[[:space:]]*$' species.ids 2>/dev/null; then
    : > "!{sample_id}.species.fa"
    echo "[ExtractSpeciesSeqs] sample=!{sample_id} species_ids=0 sequences_written=0 (no species-level)" >&2
    exit 0
  fi

  # 2) Filtrar FASTA usando normalización de headers
  awk '
    function norm_once(s) {
      gsub(/^[ \t]+|[ \t]+$/, "", s);
      sub(/^>/, "", s);
      sub(/[ \t].*$/, "", s);
      sub(/\|.*$/, "", s);
      sub(/\/.*$/, "", s);
      sub(/\.[0-9]+$/, "", s);
      return s;
    }
    function norm(s) { return norm_once(norm_once(s)); }

    BEGIN{
      while ( (getline line < "species.ids") > 0 ) {
        if (line ~ /^[ \t]*$/) continue;
        id = norm(line);
        if (id!="") ids[id]=1;
      }
      close("species.ids");
      keep=0;
    }
    /^>/{
      hdr = substr($0,2);
      keep = (norm(hdr) in ids);
    }
    { if (keep) print $0 }
  ' "!{fasta_file}" > "!{sample_id}.species.fa"

  n_ids=$(grep -cvE '^[[:space:]]*$' species.ids || true)
  n_seq=$(grep -c '^>' "!{sample_id}.species.fa" || true)
  echo "[ExtractSpeciesSeqs] sample=!{sample_id} species_ids=${n_ids} sequences_written=${n_seq}" >&2
  '''
}
