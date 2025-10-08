#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fq = "*.fastq"
params.outdir = "results"
params.noSpeciesPolishing = false
params.lowAbundanceThreshold = 0.02
params.isDemultiplexed = false
params.porechop_extra_end_trim = 0
params.noNanoplot = false
params.nanofilt_quality = 10
params.nanofilt_length = 1000
params.nanofilt_maxlength = 1700
params.nanofilt_headcrop = 0
params.nanofilt_tailcrop = 0
params.yacrd_c = 4
params.yacrd_n = 0.4
params.removeChimeras = false
params.megan_lcaAlgorithm = "naive"
params.megan_lcaCoveragePercent = 100
params.megan_topPercent = 10
params.megan_topPercentPolish = 5
params.megan_minPercentReadCover = 70
params.minimap2_k = 15
params.minimap2_f = 1000
params.minimap2_x = "map-ont"
params.minimap2_KM = 200
params.help = false

params.silvaFasta = "./silvadb/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
params.silvaTaxNcbiSp = "./silvadb/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz"
params.silvaTaxmap = "./silvadb/Exports/taxonomy/ncbi/taxmap_slv_ssu_ref_nr_138.1.txt.gz"

params.silvaFastaURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
params.silvaTaxNcbiSpURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz"
params.silvaTaxmapURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz"

params.fullSilva = false

sayHi()

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validation of parameters
def parameters_expected = [
  'help',
  'fq', 
  'outdir',
  'noSpeciesPolishing','no-species-polishing',
  'lowAbundanceThreshold','low-abundance-threshold',
  'isDemultiplexed', 'is-demultiplexed',
  'porechop_extra_end_trim',
  'noNanoplot', 'no-nanoplot',
  'nanofilt_quality',
  'nanofilt_length',
  'nanofilt_maxlength',
  'nanofilt_headcrop',
  'nanofilt_tailcrop',
  'yacrd_c',
  'yacrd_n',
  'removeChimeras', 'remove-chimeras',
  'megan_lcaAlgorithm', 'megan_lca-algorithm',
  'megan_lcaCoveragePercent', 'megan_lca-coverage-percent',
  'megan_topPercent', 'megan_top-percent',
  'megan_topPercentPolish', 'megan_top-percent-polish',
  'megan_minPercentReadCover', 'megan_min-percent-read-cover',
  'minimap2_k',
  'minimap2_f',
  'minimap2_x',
  'minimap2_KM',
  'silvaFasta', 'silva-fasta',
  'silvaTaxNcbiSp', 'silva-tax-ncbi-sp',
  'silvaTaxmap', 'silva-taxmap',
  'silvaFastaURL', 'silva-fasta-URL',
  'silvaTaxNcbiSpURL', 'silva-tax-ncbi-sp-URL',
  'silvaTaxmapURL', 'silva-taxmap-URL',
  'fullSilva', 'full-silva'
  ] as Set

def parameter_diff = params.keySet() - parameters_expected
if (parameter_diff.size() != 0){
   exit 1, "[Pipeline error] Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
}

helloParameters()

if (! params.fq ){
  println("You must provide at least 1 fastq file using --fq flag.")
  System.exit(1)
}

if (!["naive", "weighted", "longReads"].contains(params.megan_lcaAlgorithm)){
  println("lcaAlgorithm not valid, must be one of 'naive', 'weighted', or 'longReads'")
  System.exit(1)
}

Channel
  .fromPath(params.fq, checkIfExists: true)
  .ifEmpty { exit 1, "Non fastq files found: ${params.fq}" }
  .set{ fqs_ch }

// include modules
include {GetVersions} from './modules/processes'
include {Fastq2Fasta} from './modules/processes'
include {MakeDB} from './modules/processes'
include {Minimap2} from './modules/processes'
include {MeganLca} from './modules/processes'
include {GetReadInfo} from './modules/processes'
include {ComputeAbundances} from './modules/processes'
include {Polish} from './workflows/PolishMinimap'

// include sub-workflows
include {SetSilva} from './workflows/Silva'
include {Demultiplex} from './workflows/Demultiplex'
include {QFilt} from './workflows/QFiltWorkflow'
include {QCheck} from './workflows/QCheckWorkflow'

/*
 * NUEVO: Extrae todas las lecturas clasificadas a nivel especie y genera ${sample}.species.fa
 * Asume que el archivo read_info contiene columnas de ID (read_id/readid/id/read)
 * y de rango (rank/tax_rank/level), y filtra filas cuyo rango sea "species".
 */
process ExtractSpeciesSeqs {
  tag "$sample_id"

  input:
    tuple val(sample_id), path(fasta_file), path(readinfo_file)

  output:
    tuple val(sample_id), path("${sample_id}.species.fa")

  shell:
  '''
  set -euo pipefail

  # 1) Obtener IDs con rango species (autodetección de delimitador y columnas)
  python3 - << 'PY' > species.ids
import sys, csv, re

readinfo = "${readinfo_file}"
rank_headers = re.compile(r'^(rank|tax_?rank|level)$', re.IGNORECASE)
species_pat  = re.compile(r'^species$', re.IGNORECASE)

with open(readinfo, newline="") as f:
    sample = f.read(2048)
    f.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\\t;")
    except Exception:
        class D(csv.Dialect):
            delimiter="\\t"; quotechar='"'; doublequote=True; lineterminator="\\n"; skipinitialspace=True
        dialect = D()
    reader = csv.reader(f, dialect)
    rows = list(reader)

if not rows:
    sys.exit(0)

header = [h.strip() for h in rows[0]]
hl = [h.lower() for h in header]

# Columna de ID (prioriza campos habituales)
id_idx = None
for cand in ('read_id','readid','id','read'):
    if cand in hl:
        id_idx = hl.index(cand)
        break
if id_idx is None:
    id_idx = 0  # fallback

# Columna de RANK
rank_idx = None
for i,h in enumerate(header):
    if rank_headers.match(h.strip()):
        rank_idx = i
        break
if rank_idx is None:
    rank_idx = 1 if len(header) > 1 else 0

for row in rows[1:]:
    if not row:
        continue
    try:
        r = row[rank_idx].strip()
    except IndexError:
        continue
    if species_pat.match(r):
        try:
            rid = row[id_idx].strip()
        except IndexError:
            continue
        if rid:
            print(rid)
PY

  # 2) Filtrar del FASTA solo esas lecturas (sin dependencias externas)
  awk '
    BEGIN{
      while((getline line < "species.ids")>0){ids[line]=1}
      close("species.ids")
    }
    /^>/{
      split(substr($0,2), a, /[ \\t]/)
      curr = a[1]
      keep = (curr in ids)
    }
    { if(keep) print $0 }
  ' "${fasta_file}" > "${sample_id}.species.fa"

  # Asegurar existencia del archivo (aunque esté vacío)
  touch "${sample_id}.species.fa"
  '''
}

workflow {
  GetVersions()

  // Download, if not available, and generate synonyms
  SetSilva()
    SetSilva.out.fasta
      .set{ silva_fasta_ch }
    SetSilva.out.synonyms
      .set{ silva_synonyms_ch }

  // Demultiplex if not already
  if (! params.isDemultiplexed ){
    Demultiplex( fqs_ch )
    Demultiplex.out
      .flatten()
      .map { file -> tuple(file.baseName, file) }
      .set{ barcode_ch }
  } else {
    fqs_ch
      .flatten()
      .map { file -> tuple(file.baseName, file) }
      .set{ barcode_ch }
  }

  // Quality control sub-workflow
  QFilt( barcode_ch )
  QFilt.out.set{ filtered_scrubbed_ch }
  if (! params.noNanoplot ) {
    QCheck( barcode_ch, filtered_scrubbed_ch )
  }

  Fastq2Fasta( filtered_scrubbed_ch )
  Fastq2Fasta.out.set{ fasta_ch }

  MakeDB( silva_fasta_ch )
  Minimap2( fasta_ch, MakeDB.out )
  MeganLca( Minimap2.out, silva_synonyms_ch )
  MeganLca.out
    .collectFile(storeDir: "$params.outdir/Rma"){
      val, file ->
      ["${val}.rma", file]
    }
  GetReadInfo( MeganLca.out )
  GetReadInfo.out
      .set{ base_read_assingments_ch }

  // Publish read taxonomy assignments
  base_read_assingments_ch
      .collectFile(storeDir: "$params.outdir/Read_Assignments") {
          val, file -> 
          [ "${val}.read_info" , file ]
      }

  // Compute taxa counts and classification
  base_read_assingments_ch
    .map{val, file -> file}
    .collect()
    .set{ all_read_assignments }
  
  ComputeAbundances( 
    all_read_assignments, 
    silva_synonyms_ch, 
    !params.noSpeciesPolishing 
  )

  // Publish taxa counts and classification
  ComputeAbundances.out.counts
    .collectFile(storeDir: "$params.outdir/")
  ComputeAbundances.out.taxcla
    .collectFile(storeDir: "$params.outdir/")

  // NUEVO: Emparejar FASTA por muestra con read_info y extraer secuencias a nivel especie
  fasta_ch
    .join(base_read_assingments_ch)          // -> tuple(sample_id, fasta_file, readinfo_file)
    .set{ fasta_readinfo_ch }

  ExtractSpeciesSeqs( fasta_readinfo_ch )

  ExtractSpeciesSeqs.out
    .collectFile(storeDir: "$params.outdir/Species_Seqs") { sample_id, file ->
      [ "${sample_id}.species.fa", file ]
    }

  // Polish sub-Workflow
  if (!params.noSpeciesPolishing){
      Polish( 
        fasta_ch, 
        silva_fasta_ch, 
        silva_synonyms_ch, 
        base_read_assingments_ch,
        ComputeAbundances.out.silva_ids,
        ComputeAbundances.out.read_ids 
      )

      // Publish polished taxa counts and classification
      Polish.out.counts
        .collectFile(storeDir: "$params.outdir/"){
          file -> [ "COUNTS_polished.tsv" , file ]
        }
      Polish.out.taxcla
        .collectFile(storeDir: "$params.outdir/"){
          file -> [ "TAXCLA_polished.tsv" , file ]
        }
  }
}



def sayHi(){
  log.info """ ____   __  ____  ____  ____  __  __    ____ 
(  _ \\ /  \\(  _ \\(  __)(  __)(  )(  )  (  __)
 ) __/(  O ))   / ) _)  ) _)  )( / (_/\\ ) _) 
(__)   \\__/(__\\_)(____)(__)  (__)\\____/(____)
---------------------------------------------
... A full-length 16S profiling Pipeline ....
---------------------------------------------"""
}


def helpMessage() {
    log.info """
    Usage:
    A typical command for running the pipeline would be as follows:

        nextflow run microgenlab/porefile --fq "data/*.fastq"

    Input fastq file(s):
        --fq                          Path to input data (must be surrounded with quotes), e.g.: "./path/to/fastqs/*.fastq".

    Other:
        --silvaFasta                  Path to SILVA_*_SSURef_NR99_tax_silva.fasta.gz file. You can provide it 
                                      either compressed (.gz) or not. If not provided, the workflow automatically
                                      adds a download step (you must have internet connection).
        --silvaFastaURL               URL to SILVA_*_SSURef_NR99_tax_silva.fasta.gz file. It will be used if you
                                      don't provide the --silvaFasta parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz'.

        --silvaTaxNcbiSp              Path to tax_ncbi-species_ssu_ref_nr99_*.txt.gz file. You can provide it
                                      either compressed (.gz) or not. If not provided, the workflow automatically
                                      adds a download step.
        --silvaTaxNcbiSpURL           URL to tax_ncbi-species_ssu_ref_nr99_*.txt.gz file. It will be used if you
                                      don't provide the --silvaFasta parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz'.

        --silvaTaxmap                 Path to taxmap_slv_ssu_ref_nr_*.txt.gz file. You can provide it
                                      either compressed (.gz) or not. If not provided, the workflow automatically
                                      adds a download step.
        --silvaTaxmapURL              URL to taxmap_slv_ssu_ref_nr_*.txt.gz file. It will be used if you
                                      don't provide the --silvaFasta parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz'.

        --fullSilva                   By default, porefile reduces SILVA to prokatyote SSU (16S). Use this flag
                                      to deactivate the reducing step and use the full SILVA database.

        --outdir                      Name of the results directory. Default: "results".
        

    Process specific parameters:
        Porechop parameters:
        --porechop_extra_end_trim      The '--extra_end_trim' parameter of Porechop. Default: 0.

        NanoFilt parameters:
        --nanofilt_quality            The '--quality' parameter of NanoFilt. Default: 8.
        --nanofilt_length             The '--length' parameter of NanoFilt (minimum length). Default: 1000.
        --nanofilt_maxlength          The '--maxlength' parameter of NanoFilt. Default: 1700.
        --nanofilt_headcrop           The '--headcrop' parameter of NanoFilt. Default: 0.
        --nanofilt_tailcrop           The '--tailcrop' parameter of NanoFilt. Default: 0.

        Yacrd parameters:
        --yacrd_c                     The '-c' parameter of Yacrd (minimum coverage). Default: 4 .
        --yacrd_n                     The '-n' parameter of Yacrd (minimum coverage of read). Default: 0.4 .

        Minimap2 parameters:
        --minimap2_k                  The '-k' parameter of minimap2. Default: 15.
        --minimap2_x                  The '-x' parameter of minimap2. Default: 'map-ont'. Possible values: 'map-ont', 
                                      'asm5', 'asm10', 'asm20', 'map-pb', or 'map-hifi'. 
        --minimap2_f                  The '-f' parameter of minimap2. Default: 1000. Only applied in the Automap module.
        --minimap2_KM                 The '-K' parameter of minimap2, in Megabases. Default: 200.

        Megan6 parameters:
        --megan_lcaAlgorithm          The '--lcaAlgorithm' parameter of sam2rma tool (Megan6). Default: 'naive'.
                                      Possible values are: 'naive', 'weighted', or 'longReads'.
        --megan_topPercent            The '--topPercent' parameter of sam2rma tool (Megan6). Default: 10.
        --megan_topPercentPolish      The '--topPercent' parameter of sam2rma tool (Megan6) applied when polishing step
                                      is activated. Default: 5.
        --megan_minPercentReadCover   The '--minPercentReadCover' parameter of sam2rma and blast2rma tools (Megan6).
                                      Default: 70.
        --megan_lcaCoveragePercent    The '--lcaCoveragePercent' parameter of sam2rma and blast2rma tools (Megan6). 
                                      Default: 100.


    Other control options:
        --isDemultiplexed             Set this flag to avoid Demultiplex sub-workflow. If set, each fastq file is 
                                      processed as a different barcode.
        --removeChimeras              Set this flag to activate the chimera-removing step with Yacrd.
        --noNanoplot                  Set this flag to avoid QCheck sub-workflow. 
        --noSpeciesPolishing          Avoid the polishing sub-workflow.
        --lowAbundanceThreshold       The threshold of total abundance (from 0 to 1) to be considered as "low", and
                                      which the pipeline will try to re assign (default: 0.02).

    Container options (note single dash usage!):
        -profile docker               Use docker as container engine (default).
        -profile singularity          Use singularity as container engine.
        -profile podman               Use podman as container engine.

    Help:
        --help                        Print this help and exit.

    Authors: 
        Cecilia Salazar (csalazar@pasteur.edu.uy) 
        Ignacio Ferrés (igferres@gmail.com)
        Matías Giménez (mgimenez@pasteur.edu.uy)
        Gregorio Iraola (giraola@pasteur.edu.uy)

    Microbial Genomics Laboratory
    Institut Pasteur de Montevideo (Uruguay)

    """.stripIndent()
}


def helloParameters(){

  log.info """  Nextflow-version:             $nextflow.version
  Porefile-version:             $workflow.manifest.version
  Porefile-commit:              $workflow.commitId
  Porefile-branch:              $workflow.revision
  Profile:                      $workflow.profile
  Work directory:               $workflow.workDir
  Container:                    $workflow.container
  Input directory:              $params.fq
  Output directory:             $params.outdir
  Command:                      $workflow.commandLine
  ⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻""".stripIndent()
  log.info """ Minimap2 related parameters:
  --minimap2_k:                 $params.minimap2_k
  --minimap2_f:                 $params.minimap2_f
  --minimap2_x:                 $params.minimap2_x
  --minimap2_KM:                $params.minimap2_KM
⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻""".stripIndent()
  log.info """ SILVAdb related parameters: """
  if (params.silvaFasta && new File(params.silvaFasta as String).exists()){
    log.info """  --silvaFasta:                 $params.silvaFasta"""
  }else{
    log.info """ SILVAdb fasta file not provided. Download URL:
   --silvaFastaURL           $params.silvaFastaURL""".stripIndent()
  }

  if (params.silvaTaxNcbiSp && new File(params.silvaTaxNcbiSp as String).exists()){
    log.info """  --silvaTaxNcbiSp:             $params.silvaTaxNcbiSp"""
  }else{
    log.info """SILVAdb tax_ncbi-species file not provided. Download URL:
  --silvaTaxNcbiSpURL       $params.silvaTaxNcbiSpURL""".stripIndent()
  }

  if (params.silvaTaxmap && new File(params.silvaTaxmap as String).exists()){
    log.info """  --silvaTaxmap:                $params.silvaTaxmap"""
  }else{
    log.info """SILVAdb taxmap file not provided. Download URL:
  --silvaTaxmapURL          $params.silvaTaxmapURL""".stripIndent()
  }

  if (params.fullSilva){
    log.info """ Full SILVAdb selected:""".stripIndent()
  }else{
    log.info """ Reduce SILVAdb selected:""".stripIndent()
  }
  log.info """  --fullSilva:                  $params.fullSilva
⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻""".stripIndent()

  log.info """ Other process parameters:
Data is already demultiplexed?
  --isDemultiplexed:            $params.isDemultiplexed
Porechop
  --porechop_extra_end_trim:    $params.porechop_extra_end_trim
NanoFilt
  --nanofilt_quality:           $params.nanofilt_quality
  --nanofilt_length:            $params.nanofilt_length
  --nanofilt_maxlength:         $params.nanofilt_maxlength
  --nanofilt_headcrop:          $params.nanofilt_headcrop
  --nanofilt_tailcrop:          $params.nanofilt_tailcrop
Yacrd
  --removeChimeras:             $params.removeChimeras
  --yacrd_c:                    $params.yacrd_c
  --yacrd_n:                    $params.yacrd_n
NanoPlot
  --noNanoplot:                 $params.noNanoplot
⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻""".stripIndent()
}
