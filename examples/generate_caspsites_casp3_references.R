library(dagLogo)
library(dplyr)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  }
  normalizePath(getwd())
}

SCRIPT_DIR <- get_script_dir()
PACKAGE_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
DATASET_PATH <- file.path(PACKAGE_ROOT, "validation_data", "caspsites_casp3", "caspsites_casp3.csv")
detect_default_fasta_path <- function(package_root) {
  candidates <- c(
    file.path(package_root, "human_uniprot_2026_01", "UniProt_Proteome_Human_2026_01.fasta"),
    file.path(normalizePath(file.path(package_root, ".."), mustWork = FALSE), "human_uniprot_2026_01", "UniProt_Proteome_Human_2026_01.fasta")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }
  ""
}

DEFAULT_FASTA_PATH <- detect_default_fasta_path(PACKAGE_ROOT)
CANONICAL_AA_PATTERN <- "[^ACDEFGHIKLMNPQRSTVWY]"
WINDOW_COLUMNS <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
NUM_SUBSAMPLES <- 300L
RANDOM_SEED <- 404L
SUPPORTED_TEST_TYPES <- c("ztest", "fisher")

args <- commandArgs(trailingOnly = TRUE)
TEST_TYPE <- if (length(args) >= 1 && args[[1]] %in% SUPPORTED_TEST_TYPES) {
  args[[1]]
} else {
  Sys.getenv("PYDAGLOGO_TEST_TYPE", unset = "ztest")
}
FASTA_PATH <- if (length(args) >= 1 && !(args[[1]] %in% SUPPORTED_TEST_TYPES)) {
  args[[1]]
} else if (length(args) >= 2) {
  args[[2]]
} else {
  Sys.getenv("PYDAGLOGO_FASTA", unset = DEFAULT_FASTA_PATH)
}
OUTPUT_ROOT <- file.path(
  PACKAGE_ROOT,
  "validation_data",
  "caspsites_casp3",
  if (identical(TEST_TYPE, "fisher")) "reference_bundle_fisher" else "reference_bundle"
)

if (!file.exists(DATASET_PATH)) {
  stop("Public Caspase-3 dataset not found at: ", DATASET_PATH)
}
if (!TEST_TYPE %in% SUPPORTED_TEST_TYPES) {
  stop("Unsupported test type: ", TEST_TYPE)
}
if (!nzchar(FASTA_PATH) || !file.exists(FASTA_PATH)) {
  stop(
    "Provide a valid human proteome FASTA as the first argument or via PYDAGLOGO_FASTA.\n",
    "Example: Rscript examples/generate_caspsites_casp3_references.R ztest /path/to/UniProt_Proteome_Human.fasta"
  )
}

clean_residue <- function(value) {
  if (is.na(value)) {
    return("")
  }
  gsub(CANONICAL_AA_PATTERN, "", toupper(trimws(as.character(value))))
}

build_public_windows <- function(frame) {
  frame %>%
    dplyr::mutate(
      window = apply(
        as.matrix(frame[, WINDOW_COLUMNS, drop = FALSE]),
        1,
        function(row_values) {
          paste0(vapply(row_values, clean_residue, character(1)), collapse = "")
        }
      )
    ) %>%
    dplyr::filter(nchar(window) == length(WINDOW_COLUMNS))
}

make_unmapped_dagpeptides <- function(seq, upstreamOffset = 3L, downstreamOffset = 4L, type = "fasta") {
  seq <- as.character(seq)
  width <- unique(nchar(seq))
  if (length(width) != 1) {
    stop("seq must contain characters of the same length.")
  }
  center <- upstreamOffset + 1L
  anchorAA <- substring(seq, center, center)
  dat <- data.frame(
    IDs = sprintf("UNMAPPED_%04d", seq_along(seq)),
    anchorAA = anchorAA,
    anchorPos = rep(center, length(seq)),
    peptide = seq,
    anchor = anchorAA,
    upstream = substring(seq, 1L, upstreamOffset),
    downstream = substring(seq, upstreamOffset + 2L, upstreamOffset + downstreamOffset + 1L),
    stringsAsFactors = FALSE
  )
  peptide_matrix <- do.call(rbind, strsplit(seq, "", fixed = TRUE))
  methods::new(
    "dagPeptides",
    data = dat,
    peptides = peptide_matrix,
    upstreamOffset = upstreamOffset,
    downstreamOffset = downstreamOffset,
    type = type
  )
}

collect_whole_proteome_matches <- function(proteome, width) {
  pattern <- paste0("([A-Z]{", width, "})")
  matches <- gregexpr(pattern, proteome@proteome$SEQUENCE)
  unlist(regmatches(proteome@proteome$SEQUENCE, matches))
}

testDAU_fixed <- function(dagPeptides,
                          dagBackground,
                          groupingScheme = "no",
                          bgNoise = NA,
                          method = "none") {
  if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") {
    stop("dagPeptides should be an object of dagPeptides class.", call. = FALSE)
  }
  if (missing(dagBackground) || class(dagBackground) != "dagBackground") {
    stop("dagBackground should be an object of dagBackground class.", call. = FALSE)
  }
  if (!is.na(bgNoise) && (bgNoise < 0 || bgNoise > 1)) {
    stop("The background noise is a number in the range of 0 to 1", call. = FALSE)
  }

  pmethod <- match.arg(method, stats::p.adjust.methods, several.ok = FALSE)
  exp <- dagPeptides@peptides
  bg  <- dagBackground@background

  if (!groupingScheme %in% ls(envir = dagLogo:::cachedEnv)) {
    stop("Unknown grouping scheme used!")
  }

  group <- get(groupingScheme, envir = dagLogo:::cachedEnv)$group
  if (is.null(group)) {
    AA <- get(groupingScheme, envir = dagLogo:::cachedEnv)$symbol
    coln <- as.character(AA)
  } else {
    coln <- names(get(groupingScheme, envir = dagLogo:::cachedEnv)$group)
  }

  convert <- function(x, gtype) {
    d <- dim(x)
    x <- as.character(x)
    for (i in seq_along(gtype)) {
      id <- x %in% gtype[[i]]
      name <- names(gtype)[i]
      x[id] <- name
    }
    matrix(x, nrow = d[1], ncol = d[2])
  }

  groupAA <- function(dat, group) {
    if (grepl("group", group)) {
      dat <- convert(dat, get(group, envir = dagLogo:::cachedEnv)$group)
    }
    dat
  }

  bg  <- lapply(bg, function(.bg) groupAA(.bg, groupingScheme))
  exp <- groupAA(exp, groupingScheme)

  counts <- function(mat, coln) {
    apply(mat, 2, function(.ele) {
      cnt <- table(.ele)[coln]
      names(cnt) <- coln
      cnt[is.na(cnt)] <- 0
      total <- sum(cnt)
      list(freq = cnt, percent = cnt / total)
    })
  }

  bg  <- lapply(bg, counts, coln = coln)
  exp <- counts(exp, coln = coln)

  exp_freq <- do.call(cbind, lapply(exp, function(.ele) .ele$freq))
  exp_freq[is.na(exp_freq)] <- 0
  exp_percent <- do.call(cbind, lapply(exp, function(.ele) .ele$percent))
  exp_percent[is.na(exp_percent)] <- 0
  rownames(exp_freq) <- rownames(exp_percent) <- coln

  if (dagBackground@testType == "ztest") {
    per_sample_bg_percent <- lapply(bg, function(.bg) {
      do.call(cbind, lapply(.bg, function(.ele) .ele$percent))
    })
    per_position_bg_percent <- lapply(seq_len(ncol(exp_percent)), function(i) {
      do.call(cbind, lapply(per_sample_bg_percent, function(.ele) .ele[, i]))
    })

    if (!is.na(bgNoise)) {
      rdirichlet <- function(n, alpha) {
        l <- length(alpha)
        x <- matrix(stats::rgamma(l * n, alpha), ncol = n, byrow = TRUE)
        sm <- x %*% rep(1, n)
        t(x / as.vector(sm))
      }
      bg_percent <- lapply(per_position_bg_percent, function(col_pct) {
        (1 - bgNoise) * col_pct + bgNoise * rdirichlet(nrow(col_pct), rep(1, ncol(col_pct)))
      })
    } else {
      bg_percent <- per_position_bg_percent
    }

    mu_percent <- do.call(cbind, lapply(bg_percent, function(.bg) {
      apply(.bg, 1, mean, na.rm = TRUE)
    }))
    diff_percent <- exp_percent - mu_percent
    diff_percent[is.na(diff_percent)] <- 0
    n <- dagBackground@numSubsamples
    std_percent <- sqrt(mu_percent * (1 - mu_percent) / n)
    statistics <- diff_percent / std_percent
    pvalue <- 2 * stats::pnorm(-abs(statistics))
  } else {
    bg <- bg[[1]]
    bg_freq <- do.call(cbind, lapply(bg, function(.ele) {
      .ele$freq
    }))
    mu_percent <- do.call(cbind, lapply(bg, function(.ele) {
      .ele$percent
    }))
    diff_percent <- exp_percent - mu_percent
    non_bg_freq <- -sweep(x = bg_freq, MARGIN = 2, STATS = colSums(bg_freq), FUN = "-")
    non_exp_freq <- -sweep(x = exp_freq, MARGIN = 2, STATS = colSums(exp_freq), FUN = "-")
    testResults <- mapply(function(x, X, y, Y) {
      fisher <- stats::fisher.test(
        matrix(c(x, X, y, Y), nrow = 2, byrow = TRUE),
        alternative = "two.sided"
      )
      c(pvalues = fisher$p.value, statistics = fisher$estimate)
    }, exp_freq, non_exp_freq, bg_freq, non_bg_freq, SIMPLIFY = FALSE)
    testOut <- as.data.frame(do.call(rbind, testResults))
    colnames(testOut) <- c("pvalues", "statistics")
    pvalue <- matrix(testOut$pvalues, nrow = nrow(exp_percent), byrow = FALSE)
    statistics <- matrix(testOut$statistics, nrow = nrow(exp_percent), byrow = FALSE)
  }

  rownames(diff_percent) <- rownames(pvalue) <- rownames(statistics) <- coln
  colnames(diff_percent) <- colnames(statistics) <- colnames(pvalue) <- c("AA-3", "AA-2", "AA-1", "AA0", "AA1", "AA2", "AA3", "AA4")

  if (pmethod != "none") {
    pvalue[] <- do.call("cbind", lapply(as.data.frame(pvalue), function(.x) {
      stats::p.adjust(p = .x, method = pmethod)
    }))
  }

  methods::new(
    "testDAUresults",
    group = groupingScheme,
    difference = diff_percent,
    statistics = statistics,
    pvalue = pvalue,
    background = mu_percent,
    motif = exp_percent,
    testType = dagBackground@testType,
    upstreamOffset = dagPeptides@upstreamOffset,
    downstreamOffset = dagPeptides@downstreamOffset
  )
}

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
input_frame <- read.csv(DATASET_PATH, stringsAsFactors = FALSE, check.names = FALSE)
message("Loaded public dataset rows: ", nrow(input_frame))
prepared <- build_public_windows(input_frame)
message("Prepared valid 8-mer windows: ", nrow(prepared))
message("Loading proteome from: ", FASTA_PATH)
proteome <- prepareProteome(fasta = FASTA_PATH, species = "Homo sapiens")
message("Precomputing whole-proteome ", unique(nchar(prepared$window)), "-mer background matches once.")
whole_proteome_matches <- collect_whole_proteome_matches(proteome, width = unique(nchar(prepared$window)))
message("Precomputed whole-proteome matches: ", length(whole_proteome_matches))
message("Proteome loaded. Building fast unmapped dagPeptides foreground.")
seqs_4ice <- make_unmapped_dagpeptides(
  seq = prepared$window,
  upstreamOffset = 3L,
  downstreamOffset = 4L,
  type = proteome@type
)
message("Foreground formatted. Building ", TEST_TYPE, " background model from the proteome.")
if (identical(TEST_TYPE, "ztest")) {
  bg_model <- dagLogo:::buildZTestBackgroundModel(
    dagPeptides = seqs_4ice,
    matches = whole_proteome_matches,
    numSubsamples = NUM_SUBSAMPLES,
    rand.seed = RANDOM_SEED,
    replacement = FALSE,
    testType = "ztest"
  )
} else {
  bg_model <- dagLogo:::initiateBackgroundModel(
    background = list(do.call(rbind, strsplit(whole_proteome_matches, "", fixed = TRUE))),
    numSubsamples = 1L,
    testType = "fisher"
  )
}
message("Background model built. Running ", TEST_TYPE, ".")
t_res <- testDAU_fixed(dagPeptides = seqs_4ice, dagBackground = bg_model)
message("Statistics computed. Writing reference bundle.")

write.csv(data.frame(window = prepared$window), file.path(OUTPUT_ROOT, "input_windows.csv"), row.names = FALSE)
write.csv(as.data.frame(seqs_4ice@peptides), file.path(OUTPUT_ROOT, "input_windows_aligned_matrix.csv"))
write.csv(as.data.frame(t_res@difference), file.path(OUTPUT_ROOT, "difference_raw.csv"))
write.csv(as.data.frame(t_res@statistics), file.path(OUTPUT_ROOT, "statistics_raw.csv"))
write.csv(as.data.frame(t_res@pvalue), file.path(OUTPUT_ROOT, "pvalue_raw.csv"))
write.csv(as.data.frame(t_res@background), file.path(OUTPUT_ROOT, "background_mean.csv"))
write.csv(as.data.frame(t_res@motif), file.path(OUTPUT_ROOT, "motif_observed.csv"))
background_samples_long <- bind_rows(
  lapply(seq_along(bg_model@background), function(sample_id) {
    as.data.frame(bg_model@background[[sample_id]], stringsAsFactors = FALSE) %>%
      dplyr::mutate(sample_id = sample_id, .before = 1)
  })
)
background_samples_output_path <- file.path(
  OUTPUT_ROOT,
  if (identical(TEST_TYPE, "fisher")) "background_samples.csv.gz" else "background_samples.csv"
)
if (identical(TEST_TYPE, "fisher")) {
  background_con <- gzfile(background_samples_output_path, open = "wt")
  write.csv(background_samples_long, background_con, row.names = FALSE)
  close(background_con)
} else {
  write.csv(background_samples_long, background_samples_output_path, row.names = FALSE)
}
writeLines(
  c(
    paste0("test_type=", TEST_TYPE),
    paste0("fasta_path=", normalizePath(FASTA_PATH)),
    paste0("num_subsamples=", if (identical(TEST_TYPE, "fisher")) 1L else NUM_SUBSAMPLES),
    paste0("rand_seed=", RANDOM_SEED)
  ),
  file.path(OUTPUT_ROOT, "background_provenance.txt")
)

diff_mat <- t_res@difference
p_mat <- t_res@pvalue
p_adj <- matrix(
  p.adjust(as.vector(p_mat), method = "BH"),
  nrow = nrow(p_mat),
  ncol = ncol(p_mat),
  dimnames = dimnames(p_mat)
)
keep <- (diff_mat > 0) & (p_adj < 0.05)
diff_mat[!keep] <- 0
logo_matrix <- as.data.frame(t(diff_mat))
rownames(logo_matrix) <- WINDOW_COLUMNS
write.csv(logo_matrix, file.path(OUTPUT_ROOT, "logo_matrix.csv"))
write.csv(as.data.frame(p_adj), file.path(OUTPUT_ROOT, "pvalue_bh_global.csv"))

if (identical(Sys.getenv("PYDAGLOGO_WRITE_PNG", unset = "0"), "1")) {
  png(file.path(OUTPUT_ROOT, paste0("icelogo_caspsites_casp3_", TEST_TYPE, ".png")), width = 800, height = 300, res = 96)
  t_res@difference <- diff_mat
  dagLogo(
    testDAUresults = t_res,
    type = "diff",
    pvalueCutoff = 1,
    labels = WINDOW_COLUMNS,
    fontsize = 14,
    labelRelativeToAnchor = FALSE
  )
  dev.off()
}

message("Saved public Caspase-3 pydaglogo reference bundle to: ", OUTPUT_ROOT)
message(
  "Public benchmark settings (", TEST_TYPE, "): wholeProteome background from ",
  normalizePath(FASTA_PATH),
  " with numSubsamples=",
  if (identical(TEST_TYPE, "fisher")) 1L else NUM_SUBSAMPLES,
  ", rand.seed=",
  RANDOM_SEED
)
