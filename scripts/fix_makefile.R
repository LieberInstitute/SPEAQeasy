
#  Read in the Makefile used to build WiggleTools
makefilePath = paste0(getwd(), "/WiggleTools/src/Makefile")
makefile = readLines(makefilePath)

#  Find the line listing linked libraries
isLibLine = function(fileLine) {
  (nchar(fileLine) > 5) && (substr(fileLine, 1, 5) == "LIBS=")
}
libLine = match(TRUE, sapply(makefile, isLibLine))

#  Add the additional libraries to the command-line options, which seems to
#  be required for wiggletools to properly build
makefile[libLine] = paste(makefile[libLine], "-lcrypto -llzma -lbz2")

#  Rewrite the Makefile in the same location
writeLines(makefile, makefilePath)
