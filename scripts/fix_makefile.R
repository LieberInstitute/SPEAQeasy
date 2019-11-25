
#  Read in the Makefile used to build WiggleTools
makefilePath = paste0(getwd(), "/WiggleTools/src/Makefile")
makefile = readLines(makefilePath)

#  Return the index of the line beginning with the string [key]
getLineNum = function(key) {
    isLine = sapply(makefile, function(text_line) {
        (nchar(text_line) > nchar(key)) && (substr(text_line, 1, nchar(key)) == key)
    })
    match(TRUE, isLine)
}

#  Point to the gsl headers
incLine = getLineNum("INC=")
makefile[incLine] = paste(makefile[incLine], "-I../../include")

#  Add the additional required libraries to the command-line options
libLine = getLineNum("LIBS=")
makefile[libLine] = paste(makefile[libLine], "-lcrypto -llzma -lbz2")

#  Rewrite the Makefile in the same location
writeLines(makefile, makefilePath)
