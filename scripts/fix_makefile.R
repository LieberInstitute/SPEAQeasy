
#  Read in the Makefile used to build WiggleTools
makefilePath = paste0(getwd(), "/WiggleTools/src/Makefile")
makefile = readLines(makefilePath)

#  Return the index of the line beginning with the string [key]
getLineNum = function(key) {
    which(substr(makefile, 1, nchar(key)) == key)
}

#  GSL and libBigWig were installed to a local directory, which the
#  Makefile needs to be aware of
incLine = getLineNum("INC=")
makefile[incLine] = paste(makefile[incLine], "-I../../include")

libLine = getLineNum("LIB_PATHS=")
makefile[libLine] = paste(makefile[libLine], "-L../../lib")

#  Add the additional required libraries to the command-line options
libLine = getLineNum("LIBS=")
makefile[libLine] = paste(makefile[libLine], "-lcrypto -lbz2")

#  Rewrite the Makefile in the same location
writeLines(makefile, makefilePath)
