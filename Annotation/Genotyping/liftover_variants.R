/users/ssemick/BrainSwap/Shared_Coding_SNPs_within_Preivous1000G.bed

## load libraries
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(jaffelab)
library(rtracklayer)
library(GenomicRanges)

## set up liftover chain
hg19_to_hg38_chain = import.chain("hg19ToHg38.over.chain")

## load SNPs
snpMap = read.csv("/users/ssemick/BrainSwap/Shared_Coding_SNPs_Pass2_Map.csv", 
	as.is=TRUE, row.names=1)
snpMap$chromosome = paste0("chr", snpMap$chromosome)

## add type of mutation info, from illumina
info = read.csv("/users/ssemick/BrainSwap/Shared_Coding_SNPs.csv",as.is=TRUE)
snpMap$MutationType = info$Mutation.s.[match(rownames(snpMap), info$Name)]
snpMap$MutationType = ss(snpMap$MutationType, "_")

## retain only missense mutations
snpMap = snpMap[grep("(Non|Mis)sense", snpMap$MutationType),]

## convert to GRanges
snpMapGR = GRanges(snpMap$chromosome,
	IRanges(snpMap$pos, width=1, name= snpMap$snp.name),
	allele.1 = snpMap$allele.1, allele.2 = snpMap$allele.2)

## liftover
lifted = liftOver(snpMapGR, hg19_to_hg38_chain)
which(elementNROWS(lifted) == 0)
liftedSnps = unlist(lifted)

x = snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh38, names(liftedSnps),ifnotfound="drop")
liftedSnps = liftedSnps[x$RefSNP_id] # drop those that didnt match up
identical(pos(x), start(liftedSnps)) # TRUE

## write out
export(liftedSnps, con = "common_missense_SNVs_hg38.bed")


## check to dbsnp coordinates
seqlevels(liftedSnps) = gsub("chr", "ch", seqlevels(liftedSnps))
newSnps = snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh38, liftedSnps,type="equal")
