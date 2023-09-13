#   This script is intended to be run whenever SPEAQeasy is updated to use newer
#   versions of R or its packages. First, SPEAQeasy should be tested with the
#   latest R + packages (probably installed with "jhpce"). Then, on the same
#   machine, this script can be run to generate an R object saving the names
#   and versions of each package.
#
#   'install_r_packages.R' reads in this object to install precise versions of
#   each R package. The idea is to:
#       1. Have a mostly automated process for determining which R packages
#          SPEAQeasy uses (instead of manually searching many scripts whenever
#          things are updated and hardcoding the resulting list in the
#          installation script)
#       2. Improve reproducibility/ reduce variation in exact package versions
#          installed by different SPEAQeasy users

library(sessioninfo)
library(here)

#   Grab a list of packages based on "library" calls made in any R script in
#   the repo
command = paste(
    sprintf('grep -E "library\\(.*\\)" %s/*.R', here('scripts')),
    'sed -r "s/.*library\\([\'\\"]?([^\\"\'\\)]*).*/\\1/"',
    'sort -u',
    sep = " | "
)
packages = system(command, intern = TRUE)

#   Grab a dataframe of info about just those packages
pkg_info = package_info(pkgs = "installed")
pkg_info = pkg_info[match(packages, pkg_info$package),]

saveRDS(pkg_info, here('scripts', 'r_packages', 'r_package_info.rds'))