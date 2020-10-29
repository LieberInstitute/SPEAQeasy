#  The R package 'checkpoint' scans a directory of R files to determine which
#  packages to install. In practice, some dependencies which are not loaded
#  by SPEAQeasy are still required for some R scripts. This file makes
#  'library' calls to signal the full list of non-Bioconductor packages to
#  install.

library('usethis')
library('remotes')
