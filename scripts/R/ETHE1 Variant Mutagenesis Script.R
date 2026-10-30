# ==============================================================================
# ETHE1 variant mutagenesis - function library
# ==============================================================================
#
# Generates mutant ETHE1 protein sequences from transcript-level HGVS cDNA
# variants (NM_014297.5, c. numbering) and annotates / cross-checks the result.
#
# This file only DEFINES functions and constants (no side effects on source()).
# Use scripts/run_ethe1_pipeline.R for the command-line workflow, or call
# run_ethe1_pipeline() from R.
#
# Design notes
#   * Variants are applied to the FULL transcript (5'UTR + CDS + 3'UTR), and the
#     mutant protein is translated from the annotated start codon until the first
#     in-frame stop. This matters for frameshifts and stop-loss variants, whose
#     new reading frame legitimately runs into the 3'UTR.
#   * Protein-level descriptions (p.) are DERIVED by comparing the mutant protein
#     with the wild type. The sequence is the primary output; the p. string is a
#     convenience label and is cross-checked (see Codon_Check / Input_Protein_Check).
#
# Dependencies: Biostrings, dplyr, readr, stringr, tibble
# ==============================================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})


# ==============================================================================
# 1. REFERENCE CONSTANTS
# ==============================================================================

ETHE1_DEFAULTS <- list(
  gene        = "ETHE1",
  uniprot_id  = "O95571",
  transcript  = "NM_014297.5",   # MANE Select, isoform 1 (NP_055112.2)
  cds_start   = 25L,             # transcript position of c.1 (A of ATG)
  cds_end     = 789L,            # transcript position of c.765 (last base of stop codon)

  # UniProt O95571 (254 aa). Validated on load against the translated CDS.
  wt_protein = paste0(
    "MAEAVLRVARRQLSQRGGSGAPILLRQMFEPVSCTFTYLLGDRESREAVLIDPVLETAPRDAQLIKELGLR",
    "LLYAVNTHCHADHITGSGLLRSLLPGCQSVISRLSGAQADLHIEDGDSIRFGRFALETRASPGHTPGCVT",
    "FVLNDHSMAFTGDALLIRGCGRTDFQQGCAKTLYHSVHEKIFTLPGDCLIYPAHDYHGFTVSTVEEERTL",
    "NPRLTLSCEEFVKIMGNLNLPKPQQIDFAVPANMRCGVQTPTA"
  ),

  # Residue annotations carried over from the original script. They are used only
  # to flag variants that touch these positions (Hits_Fe_Ligand / Hits_Active_Site).
  # Fe ligands: His79, His135, Asp154 (two His + one Asp, as described for the
  # hETHE1 crystal structure, PDB 4CHL). The wider active-site list is user-supplied
  # and was NOT independently re-verified - edit it if your source differs.
  fe_ligand_residues   = c(79L, 135L, 154L),
  active_site_residues = c(79L, 83L, 135L, 154L, 157L, 161L, 162L, 163L, 197L)
)

AA3 <- c(
  A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", Q = "Gln", E = "Glu",
  G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys", M = "Met", F = "Phe",
  P = "Pro", S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val", "*" = "Ter"
)

# One-letter -> three-letter (vectorised; "*" -> "Ter")
aa3 <- function(x) unname(AA3[x])


# ==============================================================================
# 2. TRANSLATION
# ==============================================================================

#' Translate a nucleotide string from its first base until the first stop codon.
#' @param seq Character; begins at the start codon, may run on into the 3'UTR.
#' @return list(protein, stop_found, orf_nt). `protein` excludes the stop;
#'   `orf_nt` is the nucleotide ORF (including the stop codon if one was found,
#'   otherwise all complete codons available).
translate_orf <- function(seq) {
  stopifnot(is.character(seq), length(seq) == 1L)
  seq <- toupper(seq)
  n <- (nchar(seq) %/% 3L) * 3L
  if (n < 3L) {
    return(list(protein = "", stop_found = FALSE, orf_nt = ""))
  }

  codons <- substring(seq, seq(1L, n, by = 3L), seq(3L, n, by = 3L))
  aa     <- unname(Biostrings::GENETIC_CODE[codons])   # NA for ambiguous codons (e.g. N)
  stop_i <- match("*", aa)

  last <- if (is.na(stop_i)) length(aa) else stop_i - 1L
  if (anyNA(aa[seq_len(last)])) {
    stop("Ambiguous or invalid codon encountered before the first stop codon.")
  }

  list(
    protein    = paste(aa[seq_len(last)], collapse = ""),
    stop_found = !is.na(stop_i),
    orf_nt     = substr(seq, 1L, if (is.na(stop_i)) n else stop_i * 3L)
  )
}


# ==============================================================================
# 3. REFERENCE LOADING + VALIDATION
# ==============================================================================

#' Load the reference transcript and verify it reproduces the wild-type protein.
#' Stops with an informative message if anything is inconsistent.
load_reference <- function(fasta_path,
                           cds_start  = ETHE1_DEFAULTS$cds_start,
                           cds_end    = ETHE1_DEFAULTS$cds_end,
                           wt_protein = ETHE1_DEFAULTS$wt_protein,
                           transcript_id = ETHE1_DEFAULTS$transcript) {

  if (!file.exists(fasta_path)) stop("Transcript FASTA not found: ", fasta_path)

  tx_set <- Biostrings::readDNAStringSet(fasta_path)
  if (length(tx_set) < 1L) stop("No sequence found in FASTA file: ", fasta_path)
  if (length(tx_set) > 1L) {
    warning("FASTA has ", length(tx_set), " records; using the first (",
            names(tx_set)[1], ").")
  }

  header     <- names(tx_set)[1]
  transcript <- toupper(as.character(tx_set[[1]]))

  if (!grepl(transcript_id, header, fixed = TRUE)) {
    warning("FASTA header ('", header, "') does not mention ", transcript_id,
            ". Make sure this is the correct transcript.")
  }
  if (nchar(transcript) < cds_end) {
    stop("Transcript is shorter (", nchar(transcript), " nt) than CDS_END (", cds_end, ").")
  }

  orf <- translate_orf(substr(transcript, cds_start, nchar(transcript)))

  if (!orf$stop_found) stop("No in-frame stop codon found after CDS_START.")

  if (!identical(orf$protein, wt_protein)) {
    wt_chars  <- strsplit(wt_protein, "", fixed = TRUE)[[1]]
    obs_chars <- strsplit(orf$protein, "", fixed = TRUE)[[1]]
    n  <- min(length(wt_chars), length(obs_chars))
    d  <- which(wt_chars[seq_len(n)] != obs_chars[seq_len(n)])
    at <- if (length(d)) paste0("first difference at residue ", d[1]) else
      paste0("lengths differ (", length(obs_chars), " vs ", length(wt_chars), " aa)")
    stop("Translated CDS does not match the wild-type protein (", at, ").\n",
         "Check CDS_START/CDS_END and the transcript FASTA.")
  }

  expected_end <- cds_start + 3L * (nchar(wt_protein) + 1L) - 1L
  if (cds_end != expected_end) {
    stop("CDS_END (", cds_end, ") is inconsistent with the ORF; expected ", expected_end,
         " (CDS includes the stop codon).")
  }

  list(
    transcript    = transcript,
    header        = header,
    transcript_id = transcript_id,
    cds_start     = as.integer(cds_start),
    cds_end       = as.integer(cds_end),
    cds           = substr(transcript, cds_start, cds_end),
    cds_length    = as.integer(cds_end - cds_start + 1L),
    wt_protein    = wt_protein
  )
}


# ==============================================================================
# 4. READING / PARSING VARIANTS
# ==============================================================================

#' Read variant strings from a text file (one per line) or a CSV/TSV column.
read_variant_lines <- function(path, hgvs_column = NULL) {
  if (!file.exists(path)) stop("Variant file not found: ", path)

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "tsv")) {
    df  <- readr::read_delim(path, delim = if (ext == "csv") "," else "\t",
                             col_types = readr::cols(.default = "c"),
                             progress = FALSE, show_col_types = FALSE)
    col <- hgvs_column
    if (is.null(col)) {
      candidates <- c("hgvs_cdna", "hgvs", "cdna", "variant")
      hit <- names(df)[tolower(names(df)) %in% candidates]
      if (length(hit) == 0L) {
        stop("Could not find an HGVS column in ", basename(path),
             ". Columns: ", paste(names(df), collapse = ", "),
             ". Specify one with hgvs_column / --hgvs-column.")
      }
      col <- hit[1]
    }
    if (!col %in% names(df)) stop("Column '", col, "' not found in ", basename(path))
    lines <- df[[col]]
  } else {
    lines <- readr::read_lines(path)
  }

  lines <- stringr::str_trim(stringr::str_replace(lines, "^\ufeff", ""))
  lines <- lines[!is.na(lines) & lines != "" & !startsWith(lines, "#")]
  lines[!tolower(lines) %in% c("hgvs", "hgvs_cdna", "cdna", "variant")]
}

#' Split a raw input line into its cDNA part, optional protein part and accession.
#' Accepts e.g. "c.488G>A", "NM_014297.5:c.488G>A",
#' "NC_000019.10(NM_014297.5):c.488G>A", "c.488G>A (p.Arg163Gln)".
split_variant_line <- function(line) {
  raw       <- stringr::str_trim(line)
  accession <- stringr::str_extract(raw, "N[MR]_[0-9]+(\\.[0-9]+)?")
  cdna      <- stringr::str_extract(raw, "(?<![A-Za-z])c\\.[^\\s,;|(]+")
  protein   <- stringr::str_extract(raw, "p\\.\\(?[^\\s,;|)]+\\)?")
  protein   <- stringr::str_replace_all(protein, "[()]", "")
  list(raw = raw, accession = accession, cdna = cdna, protein = protein)
}

# --- condition helpers ---------------------------------------------------------
ethe1_condition <- function(msg, class) {
  structure(class = c(class, "error", "condition"),
            list(message = msg, call = NULL))
}
not_modelled <- function(msg) stop(ethe1_condition(msg, "ethe1_not_modelled"))
ref_mismatch <- function(msg) stop(ethe1_condition(msg, "ethe1_ref_mismatch"))

#' Parse a c. HGVS string into a mutation description.
#' Supported: SNV, del, dup, ins, delins (optionally with explicit bases).
#' Intronic/splice (c.123+1G>A) and UTR (c.-5, c.*3) variants are reported as
#' "not modelled" rather than silently mis-handled.
parse_cdna_hgvs <- function(cdna) {
  if (is.na(cdna)) stop("No c. notation found in input line.")
  h <- cdna

  if (stringr::str_detect(h, "[0-9][+-][0-9]")) {
    not_modelled("Intronic / splice-site variant (offset notation): protein effect requires splicing prediction; not modelled.")
  }
  if (stringr::str_detect(h, "^c\\.[-*]|_[-*]")) {
    not_modelled("5'UTR / 3'UTR variant (c.- / c.* notation): no change to the annotated protein sequence modelled.")
  }

  m <- stringr::str_match(h, "^c\\.([0-9]+)(?:_([0-9]+))?(.+)$")
  if (is.na(m[1, 1])) stop("Cannot parse HGVS: ", h)

  start  <- as.integer(m[1, 2])
  ranged <- !is.na(m[1, 3])
  end    <- if (ranged) as.integer(m[1, 3]) else start
  suffix <- m[1, 4]
  if (end < start) stop("Invalid range (end < start): ", h)

  rx <- function(p) stringr::regex(p, ignore_case = TRUE)
  out <- function(type, ref = "", alt = "")
    list(type = type, start = start, end = end, ref = toupper(ref), alt = toupper(alt), text = h)

  if (!ranged) {
    s <- stringr::str_match(suffix, rx("^([ACGT])>([ACGT])$"))
    if (!is.na(s[1, 1])) return(out("snv", s[1, 2], s[1, 3]))
  }
  s <- stringr::str_match(suffix, rx("^del([ACGT]*)ins([ACGT]+)$"))
  if (!is.na(s[1, 1])) return(out("delins", s[1, 2], s[1, 3]))
  s <- stringr::str_match(suffix, rx("^del([ACGT]*)$"))
  if (!is.na(s[1, 1])) return(out("deletion", s[1, 2]))
  s <- stringr::str_match(suffix, rx("^dup([ACGT]*)$"))
  if (!is.na(s[1, 1])) return(out("duplication", s[1, 2]))
  s <- stringr::str_match(suffix, rx("^ins([ACGT]+)$"))
  if (!is.na(s[1, 1])) {
    if (!ranged || end != start + 1L) {
      stop("Insertion must be written as c.X_(X+1)insSEQ: ", h)
    }
    return(out("insertion", "", s[1, 2]))
  }

  stop("Unsupported HGVS: ", h)
}


# ==============================================================================
# 5. APPLY A MUTATION TO THE TRANSCRIPT
# ==============================================================================

#' Apply a parsed mutation to the full transcript.
#' @return list(transcript = mutated transcript, ref_check = "PASS" or NA)
apply_variant <- function(ref, mut) {
  if (mut$start < 1L || mut$end > ref$cds_length) {
    stop("Position c.", mut$start, if (mut$end != mut$start) paste0("_", mut$end) else "",
         " lies outside the CDS (c.1 to c.", ref$cds_length, ").")
  }

  tx  <- ref$transcript
  off <- ref$cds_start - 1L
  s   <- mut$start + off
  e   <- mut$end + off
  seg   <- substr(tx, s, e)
  left  <- substr(tx, 1L, s - 1L)
  right <- substr(tx, e + 1L, nchar(tx))

  # Compare any reference bases written in the HGVS string with the real sequence.
  ref_check <- NA_character_
  if (nzchar(mut$ref)) {
    if (!identical(seg, mut$ref)) {
      ref_mismatch(sprintf("Reference allele mismatch at c.%s: HGVS states %s but the transcript has %s",
                           if (mut$end != mut$start) paste0(mut$start, "_", mut$end) else mut$start,
                           mut$ref, seg))
    }
    ref_check <- "PASS"
  }

  mutated <- switch(
    mut$type,
    snv         = paste0(left, mut$alt, right),
    deletion    = paste0(left, right),
    duplication = paste0(left, seg, seg, right),            # copy placed 3' of the original
    insertion   = paste0(substr(tx, 1L, s), mut$alt, substr(tx, s + 1L, nchar(tx))),
    delins      = paste0(left, mut$alt, right),
    stop("Unsupported mutation type: ", mut$type)
  )

  list(transcript = mutated, ref_check = ref_check)
}


# ==============================================================================
# 6. PROTEIN-LEVEL DESCRIPTION
# ==============================================================================

#' Derive an HGVS-style protein description by comparing mutant with wild type.
#' Sequence differences are placed at the most 3' position (HGVS 3' rule).
#' @param frameshift TRUE if the nucleotide edit changes the reading frame.
#' @return list(p, consequence, first, last): `first`/`last` = affected residues.
describe_protein_change <- function(wt, mut, stop_found, frameshift) {
  wt_aa  <- c(strsplit(wt, "", fixed = TRUE)[[1]], "*")
  mut_aa <- c(strsplit(mut, "", fixed = TRUE)[[1]], if (stop_found) "*")
  lw <- length(wt_aa)
  lm <- length(mut_aa)
  last_res <- lw - 1L

  res <- function(p, consequence, first = NA_integer_, last = NA_integer_)
    list(p = p, consequence = consequence, first = as.integer(first), last = as.integer(last))

  if (identical(wt_aa, mut_aa)) return(res("p.=", "Synonymous"))

  n  <- min(lw, lm)
  eq <- wt_aa[seq_len(n)] == mut_aa[seq_len(n)]
  i  <- if (all(eq)) n + 1L else which(!eq)[1]
  if (i > lw) stop("Internal error: mutant protein extends the wild type with an identical prefix.")

  w <- wt_aa[i]
  m <- if (i <= lm) mut_aa[i] else NA_character_
  ter_n <- function() if (stop_found) as.character(lm - i + 1L) else "?"

  # Mutant ended (no stop) before a difference appeared: frame shifted, no stop in transcript
  if (is.na(m)) {
    return(res(sprintf("p.%s%dfs", aa3(w), i),
               if (frameshift) "Frameshift" else "Indeterminate", i, last_res))
  }

  # Stop codon itself altered -> C-terminal extension
  if (w == "*") {
    return(res(sprintf("p.Ter%d%sextTer%s", i, aa3(m), ter_n()), "Stop-loss", i, i))
  }

  # Premature stop at the first changed position
  if (m == "*") {
    return(res(sprintf("p.%s%dTer", aa3(w), i),
               if (frameshift) "Frameshift" else "Nonsense", i, last_res))
  }

  # Frameshift
  if (frameshift) {
    return(res(sprintf("p.%s%d%sfsTer%s", aa3(w), i, aa3(m), ter_n()),
               "Frameshift", i, last_res))
  }

  # In-frame: trim the common suffix (never overlapping the matched prefix)
  max_s <- min(lw, lm) - (i - 1L)
  s <- 0L
  while (s < max_s && wt_aa[lw - s] == mut_aa[lm - s]) s <- s + 1L

  wm <- if (lw - s >= i) wt_aa[i:(lw - s)] else character(0)
  mm <- if (lm - s >= i) mut_aa[i:(lm - s)] else character(0)
  lenw <- length(wm)
  lenm <- length(mm)
  mm3  <- paste(aa3(mm), collapse = "")

  if (lenw == 1L && lenm == 1L) {
    return(res(sprintf("p.%s%d%s", aa3(wm), i, aa3(mm)), "Missense", i, i))
  }

  if (lenw == 0L && lenm > 0L) {
    # Insertion between residues i-1 and i; a copy of the preceding residues = duplication
    if (i - 1L >= lenm && identical(wt_aa[(i - lenm):(i - 1L)], mm)) {
      a <- i - lenm
      b <- i - 1L
      p <- if (lenm == 1L) sprintf("p.%s%ddup", aa3(wt_aa[b]), b)
           else sprintf("p.%s%d_%s%ddup", aa3(wt_aa[a]), a, aa3(wt_aa[b]), b)
      return(res(p, "In-frame duplication", a, b))
    }
    return(res(sprintf("p.%s%d_%s%dins%s", aa3(wt_aa[i - 1L]), i - 1L, aa3(wt_aa[i]), i, mm3),
               "In-frame insertion", i - 1L, i))
  }

  if (lenw > 0L && lenm == 0L) {
    j <- i + lenw - 1L
    p <- if (lenw == 1L) sprintf("p.%s%ddel", aa3(wt_aa[i]), i)
         else sprintf("p.%s%d_%s%ddel", aa3(wt_aa[i]), i, aa3(wt_aa[j]), j)
    return(res(p, "In-frame deletion", i, j))
  }

  j <- i + lenw - 1L
  p <- if (lenw == 1L) sprintf("p.%s%ddelins%s", aa3(wt_aa[i]), i, mm3)
       else sprintf("p.%s%d_%s%ddelins%s", aa3(wt_aa[i]), i, aa3(wt_aa[j]), j, mm3)
  res(p, "In-frame delins", i, j)
}


# ==============================================================================
# 7. INDEPENDENT CROSS-CHECKS
# ==============================================================================

#' Independent codon-level check for SNVs (does not use the mutant protein sequence
#' to decide the expected residue): ref codon -> alt codon -> expected amino acid,
#' compared with what the full transcript-level pipeline produced.
snv_codon_check <- function(ref, mut, mutant_protein, stop_found, consequence) {
  codon_idx <- (mut$start - 1L) %/% 3L + 1L
  k         <- (mut$start - 1L) %% 3L + 1L
  ref_codon <- substr(ref$cds, 3L * codon_idx - 2L, 3L * codon_idx)
  alt_codon <- ref_codon
  substr(alt_codon, k, k) <- mut$alt

  ref_aa <- unname(Biostrings::GENETIC_CODE[ref_codon])
  alt_aa <- unname(Biostrings::GENETIC_CODE[alt_codon])
  wt_aa  <- c(strsplit(ref$wt_protein, "", fixed = TRUE)[[1]], "*")

  if (ref_aa != wt_aa[codon_idx]) {
    return(sprintf("FAIL: reference codon %s encodes %s, expected %s", ref_codon, ref_aa, wt_aa[codon_idx]))
  }

  mut_aa <- c(strsplit(mutant_protein, "", fixed = TRUE)[[1]], if (stop_found) "*")

  ok <- if (consequence == "Start-loss") {
    alt_aa != "M"
  } else {
    codon_idx <= length(mut_aa) && identical(mut_aa[codon_idx], alt_aa)
  }

  if (ok) "PASS" else sprintf("FAIL: codon %s->%s encodes %s but pipeline gave %s",
                              ref_codon, alt_codon, alt_aa,
                              if (codon_idx <= length(mut_aa)) mut_aa[codon_idx] else "NA")
}

#' Normalise a protein HGVS string to 3-letter form ("p.R163Q" -> "p.Arg163Gln").
normalise_protein_hgvs <- function(p) {
  if (is.na(p)) return(NA_character_)
  p <- stringr::str_replace_all(p, "[()\\s]", "")
  p <- stringr::str_replace(p, "^p\\.", "")
  p <- stringr::str_replace_all(p, "\\*", "Ter")
  m <- stringr::str_match(p, "^([A-Z])([0-9]+)([A-Z]|Ter|=)$")
  if (!is.na(m[1, 1]) && m[1, 2] %in% names(AA3)) {
    third <- if (m[1, 4] %in% c("Ter", "=")) m[1, 4] else aa3(m[1, 4])
    if (!is.na(third)) p <- paste0(aa3(m[1, 2]), m[1, 3], third)
  }
  paste0("p.", p)
}

#' Compare a protein HGVS supplied with the input against the derived one.
check_input_protein <- function(input_p, derived_p) {
  if (is.na(input_p)) return(NA_character_)
  a <- normalise_protein_hgvs(input_p)
  b <- derived_p
  if (identical(a, b)) return("PASS")
  head_rx <- "^p\\.[A-Za-z]{3}[0-9]+"
  if (grepl("fs", a) && grepl("fs", b) &&
      identical(stringr::str_extract(a, head_rx), stringr::str_extract(b, head_rx))) {
    return("PASS_PARTIAL: frameshift start matches")
  }
  paste0("MISMATCH: input ", a, " vs derived ", b)
}


# ==============================================================================
# 8. ANALYSE ONE VARIANT
# ==============================================================================

evaluate_variant <- function(parts, ref, fe_ligand, active_site) {

  # --- transcript accession ---------------------------------------------------
  note <- NA_character_
  if (!is.na(parts$accession)) {
    base_in  <- sub("\\..*$", "", parts$accession)
    base_ref <- sub("\\..*$", "", ref$transcript_id)
    if (base_in != base_ref) {
      stop("Variant is on a different transcript (", parts$accession, "); expected ", ref$transcript_id, ".")
    }
    if (parts$accession != ref$transcript_id) {
      note <- sprintf("Input transcript version %s differs from reference %s; c. numbering assumed identical (allele check applied where bases are given).",
                      parts$accession, ref$transcript_id)
    }
  }

  # --- mutate + translate -----------------------------------------------------
  mut      <- parse_cdna_hgvs(parts$cdna)
  applied  <- apply_variant(ref, mut)
  mutated  <- applied$transcript
  net      <- nchar(mutated) - nchar(ref$transcript)
  frameshift <- (net %% 3L) != 0L

  start_codon <- substr(mutated, ref$cds_start, ref$cds_start + 2L)
  wt_len      <- nchar(ref$wt_protein)

  if (start_codon != "ATG") {
    # The annotated start codon is destroyed: the real translation start is unknown.
    # We do not invent a protein; the variant is flagged and excluded from FASTA export.
    desc <- list(p = "p.Met1?", consequence = "Start-loss", first = 1L, last = wt_len)
    orf  <- list(protein = NA_character_, stop_found = NA, orf_nt = NA_character_)
    status_detail <- "Start codon destroyed; initiation at a downstream ATG is not modelled."
  } else {
    orf  <- translate_orf(substr(mutated, ref$cds_start, nchar(mutated)))
    desc <- describe_protein_change(ref$wt_protein, orf$protein, orf$stop_found, frameshift)
    status_detail <- NA_character_

    if (identical(desc$p, "p.=") && mut$type == "snv") {
      codon_idx <- (mut$start - 1L) %/% 3L + 1L
      wt_full   <- c(strsplit(ref$wt_protein, "", fixed = TRUE)[[1]], "*")
      desc$p    <- sprintf("p.%s%d=", aa3(wt_full[codon_idx]), codon_idx)
    }
  }

  codon_check <- if (mut$type == "snv") {
    snv_codon_check(ref, mut, if (is.na(orf$protein)) "" else orf$protein,
                    isTRUE(orf$stop_found), desc$consequence)
  } else NA_character_

  hit <- function(sites) {
    if (is.na(desc$first)) FALSE else any(sites >= desc$first & sites <= desc$last)
  }

  list(
    HGVS_cDNA            = mut$text,
    Variant_Type         = mut$type,
    Status               = "OK",
    Status_Detail        = status_detail,
    Consequence          = desc$consequence,
    Protein_Change       = desc$p,
    Protein_Length       = if (is.na(orf$protein)) NA_integer_ else nchar(orf$protein),
    Mutant_Sequence      = orf$protein,
    Mutant_ORF_nt        = orf$orf_nt,
    Net_Length_Change_nt = as.integer(net),
    Stop_Found           = orf$stop_found,
    First_Affected_Residue = desc$first,
    Last_Affected_Residue  = desc$last,
    Hits_Fe_Ligand       = hit(fe_ligand),
    Hits_Active_Site     = hit(active_site),
    Ref_Allele_Check     = applied$ref_check,
    Codon_Check          = codon_check,
    Input_Protein_Check  = check_input_protein(parts$protein, desc$p),
    Transcript_Note      = note
  )
}

# Column template: guarantees identical columns/types for every row.
blank_result <- function() {
  list(
    HGVS_cDNA = NA_character_, Variant_Type = NA_character_, Status = NA_character_,
    Status_Detail = NA_character_, Consequence = NA_character_, Protein_Change = NA_character_,
    Protein_Length = NA_integer_, Mutant_Sequence = NA_character_, Mutant_ORF_nt = NA_character_,
    Net_Length_Change_nt = NA_integer_, Stop_Found = NA,
    First_Affected_Residue = NA_integer_, Last_Affected_Residue = NA_integer_,
    Hits_Fe_Ligand = NA, Hits_Active_Site = NA,
    Ref_Allele_Check = NA_character_, Codon_Check = NA_character_,
    Input_Protein_Check = NA_character_, Transcript_Note = NA_character_
  )
}

analyze_variant <- function(line, ref, id,
                            fe_ligand   = ETHE1_DEFAULTS$fe_ligand_residues,
                            active_site = ETHE1_DEFAULTS$active_site_residues) {

  parts <- split_variant_line(line)
  res   <- blank_result()
  res$HGVS_cDNA <- if (is.na(parts$cdna)) parts$raw else parts$cdna

  out <- tryCatch(
    evaluate_variant(parts, ref, fe_ligand, active_site),

    ethe1_not_modelled = function(e) list(Status = "NOT_MODELLED", Status_Detail = conditionMessage(e)),
    ethe1_ref_mismatch = function(e) list(Status = "ERROR", Status_Detail = conditionMessage(e),
                                          Ref_Allele_Check = paste0("FAIL: ", conditionMessage(e))),
    error              = function(e) list(Status = "ERROR", Status_Detail = conditionMessage(e))
  )
  for (nm in names(out)) res[[nm]] <- out[[nm]]

  if (res$Status == "OK") {
    res$Variant_Type <- res$Variant_Type
  } else if (is.na(res$Variant_Type)) {
    res$Variant_Type <- tryCatch(parse_cdna_hgvs(parts$cdna)$type, error = function(e) NA_character_)
  }

  res$Model_Ready <- identical(res$Status, "OK") &&
    !is.na(res$Mutant_Sequence) && nchar(res$Mutant_Sequence) > 0L

  dplyr::bind_cols(
    tibble::tibble(Variant_ID = id, Input = parts$raw,
                   Input_Protein_HGVS = parts$protein),
    tibble::as_tibble_row(res)
  )
}


# ==============================================================================
# 9. BATCH PROCESSING, QC, EXPORT
# ==============================================================================

process_variants <- function(lines, ref,
                             fe_ligand   = ETHE1_DEFAULTS$fe_ligand_residues,
                             active_site = ETHE1_DEFAULTS$active_site_residues) {
  ids <- sprintf("Variant_%03d", seq_along(lines))
  dplyr::bind_rows(lapply(seq_along(lines), function(i) {
    analyze_variant(lines[[i]], ref, ids[[i]], fe_ligand, active_site)
  }))
}

create_qc_summary <- function(tbl) {
  tab <- function(section, x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(tibble::tibble(Section = character(), Metric = character(), Count = integer()))
    t <- table(x)
    tibble::tibble(Section = section, Metric = names(t), Count = as.integer(t))
  }
  label <- function(x) sub("[:(].*$", "", x)   # "FAIL: ..." -> "FAIL"

  dup_n <- sum(duplicated(tbl$HGVS_cDNA))

  dplyr::bind_rows(
    tibble::tibble(Section = "Overview",
                   Metric = c("Total input variants", "Duplicate HGVS entries", "Model-ready protein sequences"),
                   Count  = c(nrow(tbl), dup_n, sum(tbl$Model_Ready, na.rm = TRUE))),
    tab("Status", tbl$Status),
    tab("Variant type", tbl$Variant_Type),
    tab("Consequence", tbl$Consequence),
    tab("Reference allele check", trimws(label(tbl$Ref_Allele_Check))),
    tab("SNV codon check", trimws(label(tbl$Codon_Check))),
    tab("Input protein HGVS check", trimws(label(tbl$Input_Protein_Check))),
    tab("Annotated residues", c(rep("Hits Fe ligand", sum(tbl$Hits_Fe_Ligand, na.rm = TRUE)),
                                rep("Hits active-site list", sum(tbl$Hits_Active_Site, na.rm = TRUE))))
  )
}

#' Write the mutant proteins as FASTA (e.g. for AlphaFold / ColabFold input).
#' Start-loss, failed and NOT_MODELLED variants are never exported.
#' Synonymous variants (identical to WT) are skipped unless requested.
write_variant_fasta <- function(tbl, path, ref,
                                include_wt = TRUE,
                                include_synonymous = FALSE,
                                width = 60L) {
  keep <- tbl$Model_Ready & (include_synonymous | tbl$Consequence != "Synonymous")
  keep[is.na(keep)] <- FALSE
  sub_tbl <- tbl[keep, , drop = FALSE]

  seqs  <- sub_tbl$Mutant_Sequence
  names <- sprintf("%s ETHE1 %s %s [%s]", sub_tbl$Variant_ID, sub_tbl$HGVS_cDNA,
                   sub_tbl$Protein_Change, sub_tbl$Consequence)

  if (include_wt) {
    seqs  <- c(ref$wt_protein, seqs)
    names <- c(sprintf("ETHE1_WT UniProt:%s %s wild-type", ETHE1_DEFAULTS$uniprot_id, ref$transcript_id), names)
  }

  set <- Biostrings::AAStringSet(seqs)
  names(set) <- names
  Biostrings::writeXStringSet(set, filepath = path, width = width)
  invisible(nrow(sub_tbl))
}


# ==============================================================================
# 10. END-TO-END RUNNER
# ==============================================================================

run_ethe1_pipeline <- function(transcript_fasta,
                               variants_file,
                               outdir = "results",
                               cds_start = ETHE1_DEFAULTS$cds_start,
                               cds_end   = ETHE1_DEFAULTS$cds_end,
                               hgvs_column = NULL,
                               include_synonymous_in_fasta = FALSE,
                               quiet = FALSE) {

  say <- function(...) if (!quiet) message(...)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  say("== ETHE1 variant mutagenesis pipeline ==")
  ref <- load_reference(transcript_fasta, cds_start = cds_start, cds_end = cds_end)
  say("Reference validated: CDS ", ref$cds_length, " nt, protein ", nchar(ref$wt_protein), " aa (",
      ref$transcript_id, ", UniProt ", ETHE1_DEFAULTS$uniprot_id, ")")

  lines <- read_variant_lines(variants_file, hgvs_column)
  say("Variants read: ", length(lines))

  tbl <- process_variants(lines, ref)
  qc  <- create_qc_summary(tbl)

  f_tbl   <- file.path(outdir, "ETHE1_variant_sequences.csv")
  f_qc    <- file.path(outdir, "ETHE1_QC_summary.csv")
  f_fasta <- file.path(outdir, "ETHE1_AlphaFold.fasta")
  f_info  <- file.path(outdir, "ETHE1_run_info.txt")

  readr::write_csv(tbl, f_tbl)
  readr::write_csv(qc, f_qc)
  n_fasta <- write_variant_fasta(tbl, f_fasta, ref, include_synonymous = include_synonymous_in_fasta)

  writeLines(c(
    paste("Run date:      ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("R version:     ", R.version.string),
    paste("Biostrings:    ", as.character(utils::packageVersion("Biostrings"))),
    paste("Transcript:    ", ref$transcript_id, " (FASTA header: ", ref$header, ")", sep = ""),
    paste("CDS window:    ", ref$cds_start, "-", ref$cds_end, " (transcript coordinates)", sep = ""),
    paste("Transcript md5:", unname(tools::md5sum(transcript_fasta))),
    paste("Variants file: ", basename(variants_file), " (md5 ", unname(tools::md5sum(variants_file)), ")", sep = ""),
    paste0("FASTA records:  ", n_fasta, " variant sequence(s) + WT")
  ), f_info)

  # ---- console report -----------------------------------------------------------
  bad <- tbl[tbl$Status != "OK", c("Variant_ID", "HGVS_cDNA", "Status", "Status_Detail")]
  if (nrow(bad) > 0L) {
    say("\n", nrow(bad), " variant(s) not processed:")
    if (!quiet) print(as.data.frame(bad), row.names = FALSE, right = FALSE)
  }
  mism <- tbl[!is.na(tbl$Codon_Check) & !startsWith(tbl$Codon_Check, "PASS"), ]
  mism2 <- tbl[!is.na(tbl$Input_Protein_Check) & startsWith(tbl$Input_Protein_Check, "MISMATCH"), ]
  if (nrow(mism) > 0L) say("WARNING: ", nrow(mism), " SNV(s) failed the independent codon check - inspect ", f_tbl)
  if (nrow(mism2) > 0L) say("NOTE: ", nrow(mism2), " variant(s) where the input p. notation differs from the derived one - inspect ", f_tbl)

  say("\nProcessed OK: ", sum(tbl$Status == "OK"), " / ", nrow(tbl),
      " | FASTA sequences written: ", n_fasta)
  say("Outputs in: ", normalizePath(outdir))

  invisible(list(table = tbl, qc = qc, reference = ref))
}
