#  This script exists to document a semi-automated method for producing a list
#  of R packages, read by "checkpoint::checkpoint" when locally installing R
#  and required packages. At the time of producing this script,
#  "checkpoint::checkpoint" halts with an error if any Bioconductor packages
#  are present in the directory passed to "project_dir", thus motivating
#  a workaround where a dummy directory including only non-Bioc library calls
#  is used.
#
#  Towards the bottom, code is also included to produce a list of R packages
#  for use in 'scripts/check_R_packages_JHPCE.R'.

#  Pull all library calls from R scripts used by SPEAQeasy to form a list
grep -E "library\(.*\)" ../*.R | cut -d ":" -f 2 | tr  -d " " | sort -u > r_package_list.R

#  Then manually remove lines associated with Bioconductor packages and
#  "jaffelab", and add "BiocManager"

#  Also, document how all R packages were found for use in
#  'scripts/check_R_packages_JHPCE.R', since in this script there is no need to
#  treat "ordinary" and "Bioconductor" packages differently
grep -E "library\(.*\)" ../*.R \
    | cut -d ":" -f 2 \
    | tr  -d " " \
    | sed "s/'/\"/g" \
    | sort -u \
    | cut -d "(" -f 2 \
    | cut -d ")" -f 1 \
    | grep -Ev "library|jaffelab" \
    | paste -sd ","
