#!/usr/bin/env Rscript
## Restore the recorded package versions into a repository-local library.
script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(script) == 1L)
root <- dirname(normalizePath(sub("^--file=", "", script)))
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (getRversion() != "4.5.2" || Sys.info()[["sysname"]] != "Darwin" ||
    !R.version$arch %in% c("aarch64", "arm64")) {
  stop("This exact-reference installer supports R 4.5.2 on arm64 macOS. ",
       "Other platforms require independent numerical validation.")
}
library_path <- file.path(root, "renv", "library")
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_path, .libPaths()))
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", lib = library_path)
  if (!requireNamespace("renv", quietly = TRUE)) {
    stop("Could not install renv into the repository-local library: ", library_path)
  }
}
expected <- read.csv(file.path(root, "environment_versions.csv"),
                     stringsAsFactors = FALSE)
renv::restore(project = root, library = library_path,
              packages = setdiff(expected$package, "grf"),
              lockfile = file.path(root, "renv.lock"), prompt = FALSE)

sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", shQuote(path)),
                    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(output, "status")) || length(output) != 1L) {
    stop("Could not compute SHA-256: ", path)
  }
  digest <- strsplit(trimws(output), "[[:space:]]+")[[1]][1]
  if (!grepl("^[0-9a-f]{64}$", digest)) stop("Invalid SHA-256 output.")
  digest
}

build_reference_grf <- function() {
  if (!all(nzchar(Sys.which(c("clang++", "make", "shasum"))))) {
    stop("Building reference grf requires clang++, make, and shasum.")
  }
  runtime_home <- normalizePath(R.home())
  build_dir <- file.path(library_path, ".grf-build")
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
  build_dir <- normalizePath(build_dir)
  attempt <- tempfile("attempt-", tmpdir = build_dir)
  dir.create(attempt)
  source_url <- "https://cran.r-project.org/src/contrib/Archive/grf/grf_2.4.0.tar.gz"
  source_sha <- "b53fb9b27a9d2e0e4f9da46e787c317118ace9e45e659d128af2ed8ca9176bc5"
  archive <- file.path(attempt, "grf_2.4.0.tar.gz")
  download.file(source_url, archive, method = "libcurl", mode = "wb")
  if (!identical(sha256(archive), source_sha)) {
    stop("The pinned CRAN grf source checksum did not match.")
  }
  members <- utils::untar(archive, list = TRUE)
  if (!all(startsWith(members, "grf/")) ||
      any(grepl("(^|/)\\.\\.(/|$)", members))) {
    stop("Unexpected paths in the grf source archive.")
  }
  utils::untar(archive, files = "grf/DESCRIPTION", exdir = attempt)
  desc <- read.dcf(file.path(attempt, "grf", "DESCRIPTION"))
  stopifnot(desc[1, "Package"] == "grf", desc[1, "Version"] == "2.4.0")

  # Source-build child processes must use this R, not macOS's current framework.
  home <- file.path(attempt, "r-home")
  dir.create(home)
  entries <- list.files(runtime_home, all.files = TRUE, no.. = TRUE)
  entries <- setdiff(entries, "bin")
  if (!all(file.symlink(file.path(runtime_home, entries),
                        file.path(home, entries)))) stop("Cannot prepare R build home.")
  dir.create(file.path(home, "bin"))
  entries <- setdiff(list.files(file.path(runtime_home, "bin"),
                                all.files = TRUE, no.. = TRUE), "R")
  if (!all(file.symlink(file.path(runtime_home, "bin", entries),
                        file.path(home, "bin", entries)))) stop("Cannot prepare R tools.")
  launcher <- file.path(home, "bin", "R")
  writeLines(c(
    "#!/bin/sh", "set -eu",
    'R_HOME="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"',
    'R_SHARE_DIR="$R_HOME/share"', 'R_INCLUDE_DIR="$R_HOME/include"',
    'R_DOC_DIR="$R_HOME/doc"',
    "export R_HOME R_SHARE_DIR R_INCLUDE_DIR R_DOC_DIR",
    'if [ "${1:-}" = "CMD" ]; then',
    "  shift", '  exec "$R_HOME/bin/Rcmd" "$@"', "fi",
    'exec "$R_HOME/bin/exec/R" "$@"'
  ), launcher)
  Sys.chmod(launcher, "0755")
  makevars <- file.path(attempt, "Makevars")
  writeLines(c(
    paste("LIBR =", shQuote(file.path(runtime_home, "lib", "libR.dylib"))),
    "CXXFLAGS += -ffp-contract=off", "CXX17FLAGS += -ffp-contract=off"
  ), makevars)
  staged <- file.path(attempt, "library")
  dir.create(staged)
  variables <- c(
    R_LIBS = paste(library_path, .Library, sep = .Platform$path.sep),
    R_LIBS_USER = library_path, R_LIBS_SITE = "",
    R_MAKEVARS_USER = makevars, R_MAKEVARS_SITE = "/dev/null",
    R_PROFILE = "/dev/null", R_PROFILE_USER = "/dev/null",
    R_ENVIRON = "/dev/null", R_ENVIRON_USER = "/dev/null",
    R_INSTALL_VANILLA = "true", MAKEFLAGS = "-j2"
  )
  old <- Sys.getenv(names(variables), unset = NA_character_)
  on.exit({
    Sys.unsetenv(names(old)[is.na(old)])
    do.call(Sys.setenv, as.list(old[!is.na(old)]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(variables))
  log <- file.path(attempt, "install.log")
  status <- system2(launcher, c("CMD", "INSTALL", "--no-multiarch",
                               paste0("--library=", shQuote(staged)),
                               shQuote(archive)), stdout = log, stderr = log)
  if (status != 0L) stop("Reference grf source build failed; inspect ", log)
  built <- read.dcf(file.path(staged, "grf", "DESCRIPTION"))
  if (built[1, "Version"] != "2.4.0" ||
      !startsWith(built[1, "Built"], "R 4.5.2;")) {
    stop("The source build did not use the required grf/R version.")
  }
  target <- file.path(library_path, "grf")
  backup <- file.path(attempt, "previous-grf")
  target_link <- Sys.readlink(target)
  had_previous <- file.exists(target) ||
    (!is.na(target_link) && nzchar(target_link))
  if (had_previous && !file.rename(target, backup)) {
    stop("Cannot preserve the previous local grf installation.")
  }
  if (!file.rename(file.path(staged, "grf"), target)) {
    if (had_previous) file.rename(backup, target)
    stop("Cannot place the reference grf build in the local library.")
  }
  compiler <- system2("clang++", "--version", stdout = TRUE, stderr = TRUE)
  write.dcf(data.frame(
    Package = "grf", Version = "2.4.0", Source = source_url,
    SourceSHA256 = source_sha, Built = built[1, "Built"],
    Compiler = paste(compiler, collapse = "\n"),
    CXXFlags = "-ffp-contract=off", RHome = runtime_home,
    SharedLibrarySHA256 = sha256(file.path(target, "libs", "grf.so")),
    stringsAsFactors = FALSE
  ), file.path(build_dir, "grf-provenance.dcf"))
  cat("Built grf 2.4.0 from pinned CRAN source with -ffp-contract=off.\n")
  cat("Build provenance:", file.path(build_dir, "grf-provenance.dcf"), "\n")
}

build_reference_grf()
.libPaths(c(library_path, .Library), include.site = FALSE)
records <- lapply(seq_len(nrow(expected)), function(i) {
  p <- expected$package[i]
  desc <- packageDescription(p)
  if (!identical(desc$Version, expected$version[i]) ||
      !requireNamespace(p, quietly = TRUE)) {
    stop("Restored package version/loadability check failed: ", p)
  }
  native <- list.files(file.path(find.package(p), "libs"),
                       pattern = "\\.(so|dylib)$", full.names = TRUE)
  data.frame(package = p, version = desc$Version,
             Built = if (is.null(desc$Built)) "" else desc$Built,
             provider = if (is.null(desc$Repository)) "" else desc$Repository,
             path = find.package(p),
             shared_library_sha256 = paste(vapply(native, sha256, character(1)),
                                           collapse = ";"),
             stringsAsFactors = FALSE)
})
write.csv(do.call(rbind, records),
          file.path(library_path, ".grf-build", "package_builds.csv"),
          row.names = FALSE)
cat("Restored and loaded all 73 exact package versions. ",
    "reproduce.py selects this local library automatically.\n", sep = "")
