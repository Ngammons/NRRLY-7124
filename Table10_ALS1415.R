# ================================================================
# ANALYSIS 1
# Baseline: ALS1415
# Compared strains: ALS2875, ALS2911, ALS2913, ALS2930
#
# Candidate definition:
#   1. Normalize every allele against the SAME reference FASTA.
#   2. Keep calls with QUAL > 80 in a compared strain.
#   3. Exclude an allele if an equivalent normalized call with QUAL > 80
#      is present in ALS1415.
#
# This is the reference-aware normalization used for the final dissertation
# Table 10 and the sequential Table 11 workflow. It resolves sequence-equivalent indels that can be
# represented differently in separate VCF files.
#
# Base R only. No additional R packages are required.
#
# Run in RStudio:
#   Put the five VCFs and reference FASTA in one folder, set that folder as
#   the working directory, and Source this script.
#
# Or command line:
#   Rscript Analysis_1_ALS1415_baseline_normalized.R INPUT_DIR OUTPUT_DIR REFERENCE_FASTA
# ================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(0L, 3L))) {
  stop("Supply either zero arguments, or: INPUT_DIR OUTPUT_DIR REFERENCE_FASTA")
}

input_dir <- if (length(args) == 3L) args[1] else getwd()
output_dir <- if (length(args) == 3L) args[2] else file.path(input_dir, "RESULTS_Analysis1_ALS1415")
reference_fasta <- if (length(args) == 3L) args[3] else file.path(input_dir, "GCA_016859295.1_ASM1685929v1_genomic.fna")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# The script accepts either the original Galaxy filenames or simplified names.
file_options <- list(
  ALS1415 = c("Galaxy83-[SnpEff eff_ on dataset 12 and 69](1).vcf", "ALS1415.vcf"),
  ALS2875 = c("Galaxy87-[SnpEff eff_ on dataset 12 and 70].vcf", "ALS2875.vcf"),
  ALS2911 = c("Galaxy71-[SnpEff eff_ on dataset 12 and 66](1).vcf", "ALS2911.vcf"),
  ALS2913 = c("Galaxy75-[SnpEff eff_ on dataset 12 and 67].vcf", "ALS2913.vcf"),
  ALS2930 = c("Galaxy79-[SnpEff eff_ on dataset 12 and 68](1).vcf", "ALS2930.vcf")
)

resolve_file <- function(options) {
  paths <- file.path(input_dir, options)
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop("Could not find any of: ", paste(options, collapse = ", "))
  hit[1]
}

vcf_paths <- vapply(file_options, resolve_file, character(1))

BASELINE <- "ALS1415"
COMPARED <- c("ALS2875", "ALS2911", "ALS2913", "ALS2930")
QUAL_CUTOFF <- 80

write_tsv <- function(x, filename) {
  write.table(x, file.path(output_dir, filename), sep = "\t", quote = TRUE,
              row.names = FALSE, na = "")
}

# ----------------------------
# Read the reference FASTA
# ----------------------------
read_fasta <- function(path) {
  if (!file.exists(path)) stop("Reference FASTA not found: ", path)
  x <- readLines(path, warn = FALSE)
  headers <- which(startsWith(x, ">"))
  if (!length(headers)) stop("No FASTA headers found in: ", path)

  ref <- list()
  for (i in seq_along(headers)) {
    h <- headers[i]
    end <- if (i < length(headers)) headers[i + 1L] - 1L else length(x)
    name <- strsplit(sub("^>", "", x[h]), "[[:space:]]+")[[1]][1]
    seq_lines <- x[(h + 1L):end]
    seq_lines <- seq_lines[nzchar(seq_lines)]
    ref[[name]] <- toupper(paste0(seq_lines, collapse = ""))
  }
  ref
}

reference <- read_fasta(reference_fasta)

last_base <- function(x) substr(x, nchar(x), nchar(x))
first_base <- function(x) substr(x, 1L, 1L)

# Remove redundant shared suffixes and prefixes while retaining at least
# one base in REF and ALT, as required for VCF representation.
minimal_representation <- function(pos, ref_allele, alt_allele) {
  while (nchar(ref_allele) > 1L && nchar(alt_allele) > 1L &&
         last_base(ref_allele) == last_base(alt_allele)) {
    ref_allele <- substr(ref_allele, 1L, nchar(ref_allele) - 1L)
    alt_allele <- substr(alt_allele, 1L, nchar(alt_allele) - 1L)
  }

  while (nchar(ref_allele) > 1L && nchar(alt_allele) > 1L &&
         first_base(ref_allele) == first_base(alt_allele)) {
    pos <- pos + 1L
    ref_allele <- substr(ref_allele, 2L, nchar(ref_allele))
    alt_allele <- substr(alt_allele, 2L, nchar(alt_allele))
  }

  list(pos = pos, ref = ref_allele, alt = alt_allele)
}

# Reference-aware left normalization for indels.
normalize_allele <- function(contig, pos, ref_allele, alt_allele) {
  if (is.null(reference[[contig]])) stop("Contig missing from FASTA: ", contig)

  z <- minimal_representation(pos, ref_allele, alt_allele)
  pos <- z$pos
  ref_allele <- z$ref
  alt_allele <- z$alt

  if (nchar(ref_allele) != nchar(alt_allele)) {
    while (pos > 1L && last_base(ref_allele) == last_base(alt_allele)) {
      previous_base <- substr(reference[[contig]], pos - 1L, pos - 1L)
      pos <- pos - 1L
      ref_allele <- paste0(previous_base,
                           if (nchar(ref_allele) > 1L) substr(ref_allele, 1L, nchar(ref_allele) - 1L) else "")
      alt_allele <- paste0(previous_base,
                           if (nchar(alt_allele) > 1L) substr(alt_allele, 1L, nchar(alt_allele) - 1L) else "")
    }
    z <- minimal_representation(pos, ref_allele, alt_allele)
    pos <- z$pos
    ref_allele <- z$ref
    alt_allele <- z$alt
  }

  list(pos = pos, ref = ref_allele, alt = alt_allele)
}

# ----------------------------
# Read and normalize a VCF
# ----------------------------
read_and_normalize_vcf <- function(path, strain) {
  lines <- readLines(path, warn = FALSE)
  h <- which(startsWith(lines, "#CHROM\t"))
  if (length(h) != 1L) stop("Invalid VCF header: ", path)

  col_names <- strsplit(sub("^#", "", lines[h]), "\t")[[1]]
  body <- if (h < length(lines)) lines[(h + 1L):length(lines)] else character()
  body <- body[nzchar(body) & !startsWith(body, "#")]
  if (!length(body)) stop("No variant records in: ", path)

  v <- read.delim(text = paste(body, collapse = "\n"), header = FALSE,
                  sep = "\t", quote = "", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE,
                  colClasses = "character")
  names(v) <- col_names

  required <- c("CHROM", "POS", "REF", "ALT", "QUAL", "INFO")
  if (!all(required %in% names(v))) stop("Required VCF columns missing: ", path)
  if (any(grepl(",", v$ALT, fixed = TRUE))) {
    stop("Multiallelic ALT found. Split multiallelic records before running this analysis: ", path)
  }

  v$POS <- as.integer(v$POS)
  v$QUAL <- as.numeric(v$QUAL)
  if (anyNA(v$POS) || anyNA(v$QUAL)) stop("Invalid POS or QUAL in: ", path)

  # Verify that every REF allele agrees with the supplied reference FASTA.
  for (i in seq_len(nrow(v))) {
    contig <- v$CHROM[i]
    if (is.null(reference[[contig]])) stop("Contig missing from FASTA: ", contig)
    observed_ref <- substr(reference[[contig]], v$POS[i],
                           v$POS[i] + nchar(v$REF[i]) - 1L)
    if (observed_ref != v$REF[i]) {
      stop("REF mismatch for ", strain, " at ", contig, ":", v$POS[i],
           ". VCF REF=", v$REF[i], "; FASTA REF=", observed_ref)
    }
  }

  norm <- lapply(seq_len(nrow(v)), function(i) {
    normalize_allele(v$CHROM[i], v$POS[i], v$REF[i], v$ALT[i])
  })

  v$NORM_POS <- vapply(norm, function(x) x$pos, integer(1))
  v$NORM_REF <- vapply(norm, function(x) x$ref, character(1))
  v$NORM_ALT <- vapply(norm, function(x) x$alt, character(1))
  v$key <- paste(v$CHROM, v$NORM_POS, v$NORM_REF, v$NORM_ALT, sep = "|")
  v$strain <- strain

  v$ANN <- vapply(strsplit(v$INFO, ";", fixed = TRUE), function(x) {
    y <- x[startsWith(x, "ANN=")]
    if (length(y)) substring(y[1], 5L) else ""
  }, character(1))

  # If multiple raw records collapse to the same normalized allele in one VCF,
  # retain the highest-QUAL record for comparison.
  v <- v[order(v$key, -v$QUAL), , drop = FALSE]
  v <- v[!duplicated(v$key), , drop = FALSE]
  row.names(v) <- NULL
  v
}

strains <- c(BASELINE, COMPARED)
vcfs <- setNames(lapply(strains, function(s) read_and_normalize_vcf(vcf_paths[[s]], s)), strains)

quality_summary <- do.call(rbind, lapply(strains, function(s) {
  v <- vcfs[[s]]
  data.frame(strain = s,
             normalized_records = nrow(v),
             qual_gt_80 = sum(v$QUAL > QUAL_CUTOFF),
             qual_le_80 = sum(v$QUAL <= QUAL_CUTOFF),
             stringsAsFactors = FALSE)
}))
write_tsv(quality_summary, "quality_summary.tsv")

high <- lapply(vcfs, function(v) v[v$QUAL > QUAL_CUTOFF, , drop = FALSE])
base_all <- vcfs[[BASELINE]]
base_high_keys <- high[[BASELINE]]$key

# Candidate occurrences in each compared strain.
results <- list()
for (s in COMPARED) {
  x <- high[[s]]
  x <- x[!(x$key %in% base_high_keys), , drop = FALSE]

  baseline_match <- match(x$key, base_all$key)
  x$baseline_qual_if_called <- base_all$QUAL[baseline_match]
  x$baseline_status <- ifelse(
    is.na(baseline_match),
    "NO_MATCHING_NORMALIZED_BASELINE_RECORD",
    "MATCHING_BASELINE_CALL_QUAL_LE_80"
  )
  x$comparison <- paste0(s, "_vs_", BASELINE)
  results[[s]] <- x

  write_tsv(x, paste0(s, "_vs_", BASELINE, "_candidate_occurrences.tsv"))
}

combined <- do.call(rbind, results)
row.names(combined) <- NULL
write_tsv(combined, "all_candidate_occurrences.tsv")

# One row per distinct normalized allele across all four compared strains.
candidate_keys <- unique(combined$key)
first_row <- match(candidate_keys, combined$key)
master <- combined[first_row,
                   c("CHROM", "NORM_POS", "NORM_REF", "NORM_ALT", "key", "ANN"),
                   drop = FALSE]
names(master)[1:4] <- c("CHROM", "POS", "REF", "ALT")

# Actual QUAL score in each compared strain, if QUAL > 80.
for (s in COMPARED) {
  idx <- match(master$key, high[[s]]$key)
  q <- high[[s]]$QUAL[idx]
  master[[paste0(s, "_QUAL")]] <- q
  master[[paste0(s, "_present")]] <- !is.na(q)
}

# Baseline call evidence. High-quality baseline matches should be impossible here.
bidx <- match(master$key, base_all$key)
master[[paste0(BASELINE, "_QUAL_if_called")]] <- base_all$QUAL[bidx]
master[[paste0(BASELINE, "_status")]] <- ifelse(
  is.na(bidx),
  "NO_MATCHING_NORMALIZED_RECORD",
  ifelse(base_all$QUAL[bidx] > QUAL_CUTOFF,
         "MATCHING_CALL_QUAL_GT_80",
         "MATCHING_CALL_QUAL_LE_80")
)
stopifnot(!any(master[[paste0(BASELINE, "_status")]] == "MATCHING_CALL_QUAL_GT_80"))

# Checkerboard matrix: 1 = qualifying high-quality normalized candidate; 0 = not qualifying.
checkerboard <- data.frame(CHROM = master$CHROM,
                           POS = master$POS,
                           REF = master$REF,
                           ALT = master$ALT,
                           key = master$key,
                           stringsAsFactors = FALSE)
for (s in COMPARED) {
  checkerboard[[s]] <- as.integer(master[[paste0(s, "_present")]])
}

write_tsv(master, "distinct_normalized_candidate_alleles_with_QUAL.tsv")
write_tsv(checkerboard, "checkerboard_presence_matrix.tsv")

comparison_summary <- data.frame(
  compared_strain = COMPARED,
  baseline = BASELINE,
  qualifying_candidate_occurrences = vapply(results, nrow, integer(1)),
  stringsAsFactors = FALSE
)
write_tsv(comparison_summary, "comparison_summary.tsv")

# Final verified dissertation counts.
expected_occurrences <- c(ALS2875 = 8L, ALS2911 = 2L, ALS2913 = 2L, ALS2930 = 4L)
if (any(comparison_summary$qualifying_candidate_occurrences !=
        expected_occurrences[comparison_summary$compared_strain])) {
  stop("Analysis 1 counts do not match the verified normalized dissertation results.")
}
if (nrow(master) != 14L) {
  stop("Expected 14 distinct normalized candidate alleles; observed ", nrow(master))
}

writeLines(c(
  "ANALYSIS 1",
  "Baseline: ALS1415",
  "Compared strains: ALS2875, ALS2911, ALS2913, ALS2930",
  "QUAL threshold: >80",
  "Comparison uses reference-normalized CHROM/POS/REF/ALT allele identity.",
  "Candidate = high-quality normalized allele in a compared strain with no equivalent high-quality normalized allele in ALS1415.",
  "A matching ALS1415 call with QUAL <=80 is reported as baseline evidence but does not exclude the candidate under the chosen threshold rule.",
  "No VCF record does not by itself prove the reference genotype; read-level confirmation is needed to establish absence.",
  "Verified qualifying occurrences: ALS2875=8; ALS2911=2; ALS2913=2; ALS2930=4.",
  "Verified distinct normalized candidate alleles: 14."
), file.path(output_dir, "READ_ME.txt"))

print(quality_summary)
print(comparison_summary)
cat("Distinct normalized candidate alleles:", nrow(master), "\n")
cat("Outputs saved to:", normalizePath(output_dir), "\n")
