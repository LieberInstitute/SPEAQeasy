library("styler")
library("biocthis")
library("here")

#  Style R scripts used by SPEAQeasy consistently. This script is intended to
#  be invoked whenever one or more R scripts is modified during development,
#  including both .R files and .Rmd files.
style_dir(
    path = here("scripts"),
    transformers = bioc_style(),
    recursive = FALSE
)

style_dir(
    path = here("documentation"),
    transformers = bioc_style(),
    recursive = FALSE,
    filetype = "Rmd"
)
