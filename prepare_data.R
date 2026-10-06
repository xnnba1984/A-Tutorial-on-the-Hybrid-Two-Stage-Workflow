#!/usr/bin/env Rscript
## Reconstruct the original ACTG175.csv without changing or filtering any data.
## Source: https://CRAN.R-project.org/package=speff2trial (1.0.5, GPL-2).
## Base R only: read the package's text data, without installing/executing it.
## See DATA_PREPARATION.md for provenance, licensing, and command examples.

SOURCE_VERSION <- "1.0.5"
SOURCE_ARCHIVE_MD5 <- "52ffb49246e35e40386625a007819848"
SOURCE_DATA_MD5 <- "38d443f539f20e2badf99192f5855f69"
EXPECTED_CSV_MD5 <- "04ccf3b0efa39a63f36287e6093cf038"
EXPECTED_COLUMNS <- c(
  "pidnum", "age", "wtkg", "hemo", "homo", "drugs", "karnof", "oprior",
  "z30", "zprior", "preanti", "race", "gender", "str2", "strat", "symptom",
  "treat", "offtrt", "cd40", "cd420", "cd496", "r", "cd80", "cd820",
  "cens", "days", "arms"
)
SOURCE_URLS <- c(
  "https://cran.r-project.org/src/contrib/speff2trial_1.0.5.tar.gz",
  "https://cran.r-project.org/src/contrib/Archive/speff2trial/speff2trial_1.0.5.tar.gz"
)

fail <- function(...) stop(..., call. = FALSE)

assert_md5 <- function(path, expected, label) {
  actual <- unname(tools::md5sum(path))
  if (!identical(actual, expected)) {
    fail(label, " checksum mismatch. Expected MD5 ", expected,
         "; got ", actual, ". No output has been published.")
  }
}

parse_options <- function(args) {
  opts <- list(output = NULL, reference = NULL, source_tarball = NULL)
  while (length(args)) {
    option <- args[[1L]]
    if (option %in% c("--help", "-h")) {
      cat(paste(
        "Usage: Rscript prepare_data.R [options]",
        "  --output FILE          Default: data/ACTG175.csv beside this script.",
        "  --reference FILE       Require exact row/column/value agreement first.",
        "  --source-tarball FILE  Use a previously downloaded, pinned CRAN archive.",
        "  --help                 Show this help without accessing the network.",
        "",
        "Uses base R only. Downloads speff2trial 1.0.5 over HTTPS unless an archive",
        "is supplied. Validates archive, package metadata, data, and CSV checksums.",
        "Never sorts, filters, imputes, recodes, or overwrites mismatching data.",
        "An equivalent existing output is checked and left unchanged.",
        "Generated data retain their upstream terms; do not assume the repo's",
        "code license grants permission to redistribute these data.",
        sep = "\n"
      ), "\n", sep = "")
      return(NULL)
    }
    keys <- c("--output" = "output", "--reference" = "reference",
              "--source-tarball" = "source_tarball")
    if (!(option %in% names(keys))) fail("Unknown option: ", option)
    if (length(args) < 2L || !nzchar(args[[2L]]) ||
        startsWith(args[[2L]], "--")) {
      fail("Missing file path after ", option)
    }
    key <- unname(keys[[option]])
    if (!is.null(opts[[key]])) fail("Duplicate option: ", option)
    opts[[key]] <- path.expand(args[[2L]])
    args <- args[-c(1L, 2L)]
  }
  if (is.null(opts$output)) {
    script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(script) != 1L) fail("Cannot locate script; supply --output FILE.")
    root <- dirname(normalizePath(sub("^--file=", "", script), mustWork = TRUE))
    opts$output <- file.path(root, "data", "ACTG175.csv")
  }
  opts
}

fetch_archive <- function(path) {
  failures <- character()
  for (url in SOURCE_URLS) {
    message("Downloading pinned source: ", url)
    result <- tryCatch({
      status <- suppressWarnings(utils::download.file(
        url, path, mode = "wb", quiet = TRUE
      ))
      if (status != 0L || !file.exists(path) || file.info(path)$size == 0) {
        fail("Download failed or returned an empty file.")
      }
      TRUE
    }, error = function(e) conditionMessage(e))
    if (isTRUE(result)) return(url)
    failures <- c(failures, paste(url, result, sep = ": "))
  }
  fail("Could not download the pinned CRAN source. Check network/HTTPS access,\n",
       "or use --source-tarball with the verified 1.0.5 archive.\n",
       paste(failures, collapse = "\n"))
}

validate_data <- function(dat, label) {
  if (!identical(names(dat), EXPECTED_COLUMNS) ||
      !identical(dim(dat), c(2139L, 27L))) {
    fail(label, ": expected 2,139 rows and the 27 original columns in order.")
  }
  if (!all(vapply(dat, is.numeric, logical(1))) ||
      any(vapply(dat, function(x) any(is.infinite(x)), logical(1)))) {
    fail(label, ": expected numeric columns with no infinite values.")
  }
  missing <- colSums(is.na(dat))
  expected_missing <- setNames(rep(0, length(EXPECTED_COLUMNS)), EXPECTED_COLUMNS)
  expected_missing[["cd496"]] <- 797
  if (!identical(missing, expected_missing)) {
    fail(label, ": expected 797 missing cd496 values and no other missing data.")
  }
  if (anyDuplicated(dat$pidnum)) fail(label, ": duplicate participant IDs.")
  if (!all(dat$treat == as.integer(dat$arms != 0)) ||
      !all(dat$r == as.integer(!is.na(dat$cd496)))) {
    fail(label, ": inconsistent treatment or CD4 observation indicators.")
  }
}

compare_csv <- function(path, expected, label) {
  if (!file.exists(path) || dir.exists(path)) fail(label, " is not a file: ", path)
  actual <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                    na.strings = "NA", row.names = NULL),
    error = function(e) fail("Cannot read ", label, ": ", conditionMessage(e))
  )
  validate_data(actual, label)
  if (!identical(as.numeric(actual$pidnum), as.numeric(expected$pidnum))) {
    fail(label, ": participant IDs or row order differ from the public source.")
  }
  differences <- vapply(EXPECTED_COLUMNS, function(nm) {
    x <- actual[[nm]]
    y <- expected[[nm]]
    sum(xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y))
  }, numeric(1))
  if (any(differences != 0)) {
    bad <- differences[differences != 0]
    fail(label, ": exact value comparison failed (column: differing cells): ",
         paste(paste(names(bad), bad, sep = ": "), collapse = ", "),
         ". Existing files have not been replaced.")
  }
  message(label, ": PASS (2,139 ordered rows; 27 ordered columns; ",
          "57,753 cells; zero differences; numeric tolerance = 0).")
  invisible(actual)
}

write_canonical_csv <- function(dat, path) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  utils::write.csv(dat, con, row.names = FALSE, na = "NA", eol = "\n")
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_options(args)
  if (is.null(opts)) return(invisible(NULL))
  old_options <- options(scipen = 0, OutDec = ".", timeout = max(120, getOption("timeout")))
  on.exit(options(old_options), add = TRUE)
  scratch <- tempfile("actg175_source_")
  if (!dir.create(scratch)) fail("Cannot create temporary directory.")
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

  archive <- opts$source_tarball
  if (is.null(archive)) {
    archive <- file.path(scratch, paste0("speff2trial_", SOURCE_VERSION, ".tar.gz"))
    source <- fetch_archive(archive)
  } else {
    if (!file.exists(archive) || dir.exists(archive)) fail("Source archive not found: ", archive)
    source <- normalizePath(archive, mustWork = TRUE)
  }
  assert_md5(archive, SOURCE_ARCHIVE_MD5, "CRAN source archive")

  members <- c("speff2trial/DESCRIPTION", "speff2trial/MD5",
               "speff2trial/data/ACTG175.txt")
  listed <- utils::untar(archive, list = TRUE, tar = "internal")
  if (!all(members %in% listed)) fail("Required source archive members are missing.")
  status <- utils::untar(archive, files = members, exdir = scratch, tar = "internal")
  if (status != 0L) fail("Could not extract the source package data.")
  package_dir <- file.path(scratch, "speff2trial")
  meta <- read.dcf(file.path(package_dir, "DESCRIPTION"))
  expected_meta <- c(Package = "speff2trial", Version = SOURCE_VERSION,
                     License = "GPL-2", Repository = "CRAN")
  if (nrow(meta) != 1L || !all(names(expected_meta) %in% colnames(meta)) ||
      !identical(unname(meta[1L, names(expected_meta)]), unname(expected_meta))) {
    fail("Source package name, version, license, or repository does not match the pin.")
  }
  manifest <- utils::read.table(file.path(package_dir, "MD5"),
                                colClasses = "character", quote = "", comment.char = "")
  entry <- manifest[manifest[[2L]] == "*data/ACTG175.txt", 1L]
  if (!identical(entry, SOURCE_DATA_MD5)) fail("Package MD5 manifest does not match the data pin.")
  data_path <- file.path(package_dir, "data", "ACTG175.txt")
  assert_md5(data_path, SOURCE_DATA_MD5, "ACTG175 source data")

  dat <- utils::read.table(data_path, header = TRUE, check.names = FALSE,
                           stringsAsFactors = FALSE, na.strings = "NA",
                           comment.char = "", dec = ".")
  validate_data(dat, "Public source")
  candidate <- file.path(scratch, "ACTG175.csv")
  write_canonical_csv(dat, candidate)
  assert_md5(candidate, EXPECTED_CSV_MD5, "Reconstructed CSV")
  compare_csv(candidate, dat, "CSV round-trip")
  if (!is.null(opts$reference)) compare_csv(opts$reference, dat, "Reference CSV")

  if (file.exists(opts$output)) {
    compare_csv(opts$output, dat, "Existing output (no overwrite)")
    message("Existing output retained: ", normalizePath(opts$output, mustWork = TRUE))
  } else {
    parent <- dirname(opts$output)
    if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
      fail("Cannot create output directory: ", parent)
    }
    staged <- tempfile(".ACTG175_", tmpdir = parent, fileext = ".csv")
    on.exit(unlink(staged), add = TRUE)
    if (!file.copy(candidate, staged, overwrite = FALSE)) fail("Cannot stage output CSV.")
    assert_md5(staged, EXPECTED_CSV_MD5, "Staged CSV")
    if (file.exists(opts$output)) fail("Output appeared during preparation; rerun to verify it.")
    if (!file.rename(staged, opts$output)) fail("Cannot publish verified output CSV.")
    message("Created: ", normalizePath(opts$output, mustWork = TRUE))
  }
  message("Source: speff2trial ", SOURCE_VERSION, " (GPL-2); ", source)
  message("Output MD5: ", unname(tools::md5sum(opts$output)))
  message("Upstream data terms: https://CRAN.R-project.org/package=speff2trial")
  message("Preparation complete. No scientific analyses were run.")
  invisible(opts$output)
}

if (sys.nframe() == 0L) {
  tryCatch(main(), error = function(e) {
    message("prepare_data.R: ERROR: ", conditionMessage(e))
    quit(save = "no", status = 1L)
  })
}
