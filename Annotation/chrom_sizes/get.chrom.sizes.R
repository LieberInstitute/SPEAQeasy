
library(Biostrings)
ss = function(x, pattern, slot = 1, ...) {sapply(strsplit(x = x, split = pattern, ...), "[", slot)}


##################
##### hg38 #######
fa = readDNAStringSet("GRCh38.primary_assembly.genome.fa")
hg38 = data.frame(chr=names(fa), len=width(fa))
hg38$chr = as.character(hg38$chr)
hg38$chr = ss(hg38$chr, " ")
write.table(hg38, file="hg38.chrom.sizes.gencode", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


##################
##### hg19 #######
fa = readDNAStringSet("GRCh37.primary_assembly.genome.fa")
hg19 = data.frame(chr=names(fa), len=width(fa))
hg19$chr = as.character(hg19$chr)
hg19$chr = ss(hg19$chr, " ")
write.table(hg19, file="hg19.chrom.sizes.gencode", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


##################
##### mm10 #######
fa = readDNAStringSet("GRCm38.primary_assembly.genome.fa")
mm10 = data.frame(chr=names(fa), len=width(fa))
mm10$chr = as.character(mm10$chr)
mm10$chr = ss(mm10$chr, " ")
write.table(mm10, file="mm10.chrom.sizes.gencode", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


##################
##### rn6 ########
fa = readDNAStringSet("Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa")
rn6 = data.frame(chr=names(fa), len=width(fa))
rn6$chr = as.character(rn6$chr)
rn6$chr = ss(rn6$chr, " ")
write.table(rn6, file="rn6.chrom.sizes.ensembl", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


