packages <- c("grf", "ggplot2", "pROC", "PRROC", "dplyr", "data.table",
              "survival", "sandwich", "lmtest", "car", "patchwork", "scales",
              "gridExtra")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
                          ". Install the documented dependencies before running.")
dir.create("result", showWarnings = FALSE)
write.csv(data.frame(package = packages, version = vapply(packages, function(p) {
  packageDescription(p, fields = "Version")
}, character(1))), "result/package_versions.csv", row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "", "RNGkind():",
             capture.output(RNGkind())), "result/session_info.txt")
expected <- read.csv("environment_versions.csv", stringsAsFactors = FALSE)
installed <- vapply(expected$package, function(p) {
  if (requireNamespace(p, quietly = TRUE)) packageDescription(p, fields = "Version") else NA_character_
}, character(1))
different <- is.na(installed) | installed != expected$version
if (any(different)) {
  warning("Versions differ from the verified environment: ",
          paste(expected$package[different], collapse = ", "),
          ". Restore renv.lock for the recorded versions and verify new results.")
}
print(read.csv("result/package_versions.csv"))
