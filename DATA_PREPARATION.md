# ACTG 175 Data Preparation

## Run

```sh
Rscript prepare_data.R
```

This creates `data/ACTG175.csv` beside `prepare_data.R`, including when the script
is invoked by absolute path from another directory. It requires only the standard
R installation and HTTPS access to CRAN. It does not install packages, load package
code, or run analyses. All 2,139 rows and 27 columns are retained in source order;
there is no sorting, filtering, imputation, recoding, or added row-number column.

An existing CSV is compared and left unchanged if equivalent. A mismatch stops
with an error; inspect the file before moving it or choosing another output path.

```sh
Rscript prepare_data.R --output /path/to/ACTG175.csv
Rscript prepare_data.R --reference /path/to/original/ACTG175.csv
Rscript prepare_data.R --help
```

The optional reference comparison uses the exact column order, participant order,
missing-value positions, and numeric values with zero tolerance. Integer versus
double storage is immaterial. Comparison failure occurs before output creation.
New outputs use LF line endings and a canonical CSV format. Equivalent existing
files with different text formatting are retained without rewriting them.

For an offline run, first obtain the pinned source archive from the source link
below, then supply it explicitly:

```sh
Rscript prepare_data.R --source-tarball /path/to/speff2trial_1.0.5.tar.gz
```

## Source and Checks

The source is [speff2trial 1.0.5 on CRAN](https://CRAN.R-project.org/package=speff2trial),
published May 31, 2022. The script retrieves its
[official source archive](https://cran.r-project.org/src/contrib/speff2trial_1.0.5.tar.gz)
and reads `speff2trial/data/ACTG175.txt`. If the current-release URL becomes
unavailable, it tries the same version under CRAN's `Archive/speff2trial/` path.
Version 1.0.5 was not in that archive directory on October 6, 2026; the fallback
is intended for a future CRAN move, not a currently verified second copy.

The archive checksum, `DESCRIPTION` package/version/license/repository fields,
data checksum and package MD5 manifest, dimensions, column order, unique IDs,
missingness, and treatment/observation indicators are checked before export. The
canonical CSV checksum then checks every value and both row and column order.
MD5 checks are portable accidental-change checks, not digital signatures. HTTPS
is required for the automatic download. Any source change requires a new audit;
do not bypass a failed checksum or silently use a different package release.

| Artifact | MD5 |
| --- | --- |
| `speff2trial_1.0.5.tar.gz` | `52ffb49246e35e40386625a007819848` |
| `data/ACTG175.txt` within that archive | `38d443f539f20e2badf99192f5855f69` |
| Reconstructed/original CSV | `04ccf3b0efa39a63f36287e6093cf038` |

The original and reconstructed CSV files were also byte-identical in the R 4.5.2
audit on October 6, 2026 (171,275 bytes), with SHA-256
`debab27bbc6f741724b79e49a8e13e1feb8db57850b17010b0635202bfc7141a`.
Only data preparation was audited here, not downstream scientific results.

## Data Terms and Interpretation

CRAN and the source `DESCRIPTION` declare GPL-2 for the package. Its
[ACTG175 documentation](https://cran.r-project.org/web/packages/speff2trial/refman/speff2trial.html#ACTG175)
describes the dataset and cites Hammer et al. (1996), *New England Journal of
Medicine*, 335:1081-1090. No separate data-specific license or permission letter
was found in the source archive. The repository's MIT code license is not evidence
that these third-party data may be redistributed under MIT. Consult the
[upstream GPL-2 terms](https://www.r-project.org/Licenses/GPL-2) before distributing
copies; this preparation method distributes the recipe, not the patient-level data.
Keep generated CSVs and downloaded source archives out of the published code package
unless the upstream terms and required notices have been addressed separately.

`cens = 1` marks an observed event, and `r = 1` marks an observed 96-week CD4
measurement. They are different indicators. The 797 missing `cd496` values are
preserved. Treatment arm exclusions and endpoint construction remain in the
scientific analysis scripts; data preparation does not perform those steps.
