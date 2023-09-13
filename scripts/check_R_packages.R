#   These packages are necessary for the installation of other packages, and
#   thus will be installed manually
meta_packages = c("BiocManager", "remotes", "sessioninfo")
for (package in meta_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package, repos = "http://cran.us.r-project.org")
    }
}

library(sessioninfo)

#   Read in package info
pkg_info = readRDS('scripts/r_packages/r_package_info.rds')

#   Install the Bioconductor packages-- versioning of packages is fixed by
#   Bioconductor version, which is fixed by R version
BiocManager::install(
    pkg_info[pkg_info$source == "Bioconductor", 'package'], update = FALSE
)

#   Install exact versions of CRAN packages
cran_pkg_df = pkg_info[grep('CRAN', pkg_info$source),]
for (i in 1:nrow(cran_pkg_df)) {
    if (!requireNamespace(cran_pkg_df$package[i], quietly = TRUE)) {
        remotes::install_version(
            cran_pkg_df$package[i], version = cran_pkg_df$version[i],
            repos = "http://cran.us.r-project.org", update = FALSE
        )
    }
}

#   Install exact commit in version-control history for GitHub-based packages
gh_packages = sub(
    '.*\\((.*)\\).*',
    '\\1',
    pkg_info[grep('^Github', pkg_info$source), 'source'], "\\("
)
remotes::install_github(gh_packages)