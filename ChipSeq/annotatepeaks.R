#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# annotatepeaks.R
#
# Description: Annotates genomic peak BED files (e.g., from CUT&RUN / ChIP-seq)
#              with the nearest genes using ChIPseeker and TxDb annotation packages.
# Version:     1.1.0
# Author:      Eduard Perez / Epigenomics
# Date:        2025-04-24
#
# Usage:       Rscript annotate_peaks.R -i peak1.bed,peak2.bed -g mm10 -o anno1.tsv,anno2.tsv
#
# Options:
#   -i, --input     Input peak file(s) in BED format
#   -o, --output    Output file path(s) for annotated peaks..
#                   If omitted, output files are generated next to input (_annotated.tsv/csv).
#   -g, --genome    Reference genome assembly tag (e.g., rn6, mm10, hg38) [default: "rn6"].
#   -t, --tss       TSS region definition as 'upstream,downstream' relative to TSS
#                   (e.g., "-3000,3000") [default: "-3000,3000"].
#   -s, --same-strand Require peaks and annotated genes to be on the same strand [default: FALSE].
#   -f, --format    Output file format: "tsv" or "csv" [default: "tsv"].
#   -h, --help      Show this help message and exit.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse) # Command‑line option parsing
  library(ChIPseeker) # Peak annotation utilities
})

# --------------------------- Command‑line options ---------------------------
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Input BED file(s), comma‑separated (REQUIRED)"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    help = "Output file(s), comma‑separated; default = <input>_annotated"
  ),
  make_option(
    c("-g", "--genome"),
    type = "character",
    default = "rn6",
    help = "Genome tag: rn6 | mm10 | hg38  [default %default]"
  ),
  make_option(
    c("-t", "--tss"),
    type = "character",
    default = "-3000,3000",
    help = "TSS window as 'up,down' in bp   [default %default]"
  ),
  make_option(
    c("-s", "--same-strand"),
    type = "logical",
    default = FALSE,
    help = "TRUE = annotate only genes on the same strand"
  ),
  make_option(
    c("-f", "--format"),
    type = "character",
    default = "tsv",
    help = "Output format: tsv | csv         [default %default]"
  )
)

# Parse the options supplied by the user
opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------ Input checks --------------------------------
if (is.null(opt$input)) stop("Missing -i / --input option.", call. = FALSE)

# Split potential comma lists and remove surrounding whitespace
inputs <- trimws(strsplit(opt$input, ",")[[1]])
outputs <- if (!is.null(opt$output)) trimws(strsplit(opt$output, ",")[[1]]) else
  character(0)

if (length(outputs) > length(inputs))
  stop("More outputs than inputs.", call. = FALSE)

if (!opt$format %in% c("tsv", "csv"))
  stop("--format must be 'tsv' or 'csv'.", call. = FALSE)

# Convert the TSS window string to an integer vector
tss <- as.integer(trimws(strsplit(opt$tss, ",")[[1]]))
if (length(tss) != 2 || any(is.na(tss)))
  stop("TSS window must be two integers, e.g. -3000,3000.", call. = FALSE)

# ------------------------- Genome package mapping --------------------------
genome_map <- list(
  rn6 = list(txdb = "TxDb.Rnorvegicus.UCSC.rn6.refGene", org = "org.Rn.eg.db"),
  mm10 = list(
    txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    org = "org.Mm.eg.db"
  ),
  hg38 = list(txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene", org = "org.Hs.eg.db")
)
if (!opt$genome %in% names(genome_map))
  stop("Unknown genome tag: ", opt$genome, call. = FALSE)

# ---------------------- Ensure required packages exist ----------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Load (or install) each required Bioconductor package
load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

lapply(genome_map[[opt$genome]], load_pkg)

txdb <- get(genome_map[[opt$genome]]$txdb) # Transcript database
annoDB <- genome_map[[opt$genome]]$org # Org‑db for gene symbols

# ---------------------------- Peak annotation ------------------------------
for (i in seq_along(inputs)) {
  bed <- inputs[i]
  if (!file.exists(bed)) stop("File not found: ", bed, call. = FALSE)

  # Decide output filename: user‑provided or auto‑generated
  out <- if (i <= length(outputs)) {
    outputs[i]
  } else {
    ext <- ifelse(opt$format == "csv", "_annotated.csv", "_annotated.tsv")
    sub("\\.bed$", ext, bed, ignore.case = TRUE)
  }

  message(sprintf("[%d/%d] %s -> %s", i, length(inputs), bed, out))

  # Annotate peaks with ChIPseeker and convert result to data‑frame
  peak_df <- as.data.frame(
    annotatePeak(
      peak = bed,
      TxDb = txdb,
      tssRegion = tss,
      annoDb = annoDB,
      sameStrand = isTRUE(opt$same.strand)
    )
  )

  # Write the annotated peaks table
  if (opt$format == "csv") write.csv(peak_df, out, row.names = FALSE) else
    write.table(peak_df, out, sep = "\t", quote = FALSE, row.names = FALSE)
}

message("Finished ", length(inputs), " file(s)")
