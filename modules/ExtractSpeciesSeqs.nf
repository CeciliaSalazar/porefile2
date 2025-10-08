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

  # 1) Extract species-level read IDs from .read_info
  #    - Accept TAB, comma, or semicolon as separators
  #    - First column = read id
  #    - Line contains [S] anywhere
  #    - Skip header (NR==1)
  awk -F '[\\t,;]' '
    BEGIN{ OFS="\\t" }
    NR==1 { next }
    {
      line=$0
      gsub(/\r/,"",line)   # strip CR if present
      if (index(line,"[S]")) {
        id=$1
        sub(/^>/,"",id)
        gsub(/^[ \\t]+|[ \\t]+$/,"",id)
        if (id != "") print id
      }
    }
  ' "!{readinfo_file}" > species.ids

  # Save a copy of the raw ID list for debugging/inspection
  cp species.ids "!{sample_id}.species.ids"

  # If there are no IDs, emit empty fasta and exit cleanly
  if ! grep -qve '^[[:space:]]*$' species.ids 2>/dev/null; then
    : > "!{sample_id}.species.fa"
    echo "[ExtractSpeciesSeqs] sample=!{sample_id} species_ids=0 sequences_written=0 (no [S] in read_info)" >&2
    exit 0
  fi

  # 2) Filter FASTA using normalized matching
  #    Normalization rules for both ID list and FASTA headers:
  #      - trim spaces
  #      - drop leading '>'
  #      - cut at first space, then at first '|' and at first '/'
  #      - drop trailing .<digits> (e.g. .1)
  awk '
    function norm_once(s) {
      gsub(/^[ \\t]+|[ \\t]+$/, "", s);
      sub(/^>/, "", s);
      sub(/[ \\t].*$/, "", s);
      sub(/\\|.*$/, "", s);
      sub(/\\/.*$/, "", s);
      sub(/\\.[0-9]+$/, "", s);
      return s;
    }
    function norm(s) {
      # apply twice to be safe on odd headers
      return norm_once(norm_once(s));
    }
    BEGIN{
      n_ids=0;
      while ( (getline line < "species.ids") > 0 ) {
        if (line ~ /^[ \\t]*$/) continue;
        id = norm(line);
        if (id!="") { ids[id]=1; n_ids++; }
      }
      close("species.ids");
      keep=0;
    }
    /^>/{
      hdr = substr($0,2);
      h = norm(hdr);
      keep = (h in ids);
    }
    { if (keep) print $0 }
  ' "!{fasta_file}" > "!{sample_id}.species.fa"

  n_ids=$(grep -cvE '^[[:space:]]*$' species.ids || true)
  n_seq=$(grep -c '^>' "!{sample_id}.species.fa" || true)
  echo "[ExtractSpeciesSeqs] sample=!{sample_id} species_ids=${n_ids} sequences_written=${n_seq}" >&2
  '''
}
