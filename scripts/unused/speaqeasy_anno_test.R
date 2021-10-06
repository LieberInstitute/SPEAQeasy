#  Here I explore how to replace biomaRt with org.*.eg.db packages, for the
#  purpose of getting EntrezIDs for each gene associated with the counts found
#  in a given experiment (by running SPEAQeasy)

#  This example will involve human reference and hg38 gencode v32 GTF, and is
#  done only for genes (though code will be analagous for exons, junctions)

#  Modeled after relevant code segments in create_count_objects-human.R

library('org.Hs.eg.db')
library('rtracklayer')
library('jaffelab')
library('clusterProfiler')
library('here')
library('getopt')

#  After running SPEAQeasy on paired, forward-stranded 'hg38' test samples,
#  pass the working directory for the CountObjects process here. Same, for
#  'rn6' instead of 'hg38'
spec <- matrix(
    c(
        "work_dir_human", "h", 1, "logical", "Work dir for human",
        "work_dir_rat", "r", 1, "logical", "Work dir for human"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

#  Directly reference the GTF
gencodeGTF = import(here("Annotation", "RSeQC", "hg38", "gtf", "gencode.v32.annotation.gtf"), format='gtf')

#  Use a particular SPEAQeasy run just for testing
setwd(opt$work_dir_human)
manifest = read.table(
    'samples_complete.manifest', sep = ' ', header = FALSE,
    stringsAsFactors = FALSE
)
metrics = data.frame(
    'SAMPLE_ID' = manifest[, ncol(manifest)-1], stringsAsFactors = FALSE
)
opt = list(
    'organism' = 'hg38', 'experiment'='Jlab_experiment', 'prefix'='',
    'paired'=TRUE, 'stranded'='forward', 'ercc'=FALSE
)

###############################################################################
#  Code directly from SPEAQeasy
###############################################################################

gencodeGTF = import(con=list.files(pattern=".*\\.gtf"), format="gtf")
gencodeGENES = mcols(gencodeGTF)[which(gencodeGTF$type=="gene"),c("gene_id","type","gene_type","gene_name")]
rownames(gencodeGENES) = gencodeGENES$gene_id

gencodeEXONS = as.data.frame(gencodeGTF)[which(gencodeGTF$type=="exon"),c("seqnames","start","end","gene_id","exon_id")]
names(gencodeEXONS) = c("Chr","Start","End","gene_id","exon_gencodeID")
### exons in PAR regions
par_y = grep("PAR_Y",gencodeEXONS$gene_id)
gencodeEXONS$exon_gencodeID[par_y] = paste0(gencodeEXONS$exon_gencodeID[par_y],"_PAR_Y")
gencodeEXONS = gencodeEXONS[,-4]

###############
### gene counts

#  Get filenames for gene counts in sample order present in manifest
geneFn = rep('', length(metrics$SAMPLE_ID))
for (i in 1:length(metrics$SAMPLE_ID)) {
    geneFn[i] = list.files(pattern=paste0('^',metrics$SAMPLE_ID[i], '_',
                                          opt$organism, '_.*_Genes\\.counts$'))
}
names(geneFn) = metrics$SAMPLE_ID

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]

## organize gene map
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
rownames(geneMap) = geneMap$Geneid
geneMap$gencodeID = geneMap$Geneid
geneMap$ensemblID = ss(geneMap$Geneid, "\\.")
geneMap$Geneid = NULL
geneMap$gene_type = gencodeGENES[geneMap$gencodeID,"gene_type"]
geneMap$Symbol = gencodeGENES[geneMap$gencodeID,"gene_name"]

#  Get filenames for exon counts in sample order present in manifest
exonFn = rep('', length(metrics$SAMPLE_ID))
for (i in 1:length(metrics$SAMPLE_ID)) {
    exonFn[i] = list.files(pattern=paste0('^',metrics$SAMPLE_ID[i], '_',
                                          opt$organism, '_.*_Exons\\.counts$'))
}
names(exonFn) = metrics$SAMPLE_ID

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$gencodeID = exonMap$Geneid
exonMap$ensemblID = ss(exonMap$Geneid, "\\.")
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$Geneid = NULL
exonMap$gene_type = gencodeGENES[exonMap$gencodeID,"gene_type"]
exonMap$Symbol = gencodeGENES[exonMap$gencodeID,"gene_name"]

exonMap = join(exonMap, gencodeEXONS, type="left", match="first")

###############################################################################
#  Get EntrezID
###############################################################################

geneMap$EntrezID = mapIds(org.Hs.eg.db, geneMap$Symbol, "ENTREZID", "SYMBOL")

#  TODO: select the first entrezID for each symbol

###############################################################################
#  Verify GO results agree whenever a gene matches to more than one Entrez ID
###############################################################################

a = unique(gencodeGENES$gene_name)
ent_ids = mapIds(org.Hs.eg.db, a, 'ENTREZID', "SYMBOL", multiVals='CharacterList')

#  We find that 4 genes have exactly 2 Entrez IDs per Gencode ID ("symbol")
ambig_ids = ent_ids[sapply(ent_ids, length) > 1]

#  Get character vectors of each Entrez ID set
first_ids = sapply(ambig_ids, function(x) x[1])
second_ids = sapply(ambig_ids, function(x) x[2])

#  Conveniently, elements of 'a' are NA whenever elements in 'b' are not
#  (GO results never disagree-- in a trivial sense-- and only complement 
#  eachother)
a = mapIds(org.Hs.eg.db, first_ids, 'GO', 'ENTREZID', multiVals='CharacterList')
b = mapIds(org.Hs.eg.db, second_ids, 'GO', 'ENTREZID', multiVals='CharacterList')

######################################################
#  For rat, things are more complicated
######################################################

setwd(opt$work_dir_rat)

manifest = read.table('samples_complete.manifest', sep = ' ',
    header = FALSE, stringsAsFactors = FALSE)
metrics = data.frame('SAMPLE_ID' = manifest[, ncol(manifest)-1],
    stringsAsFactors = FALSE)
opt = list('organism' = 'rn6', 'experiment'='Jlab_experiment', 'prefix'='',
    'paired'=TRUE, 'stranded'='forward', 'ercc'=FALSE)

geneFn = rep('', length(metrics$SAMPLE_ID))
for (i in 1:length(metrics$SAMPLE_ID)) {
    geneFn[i] = list.files(pattern=paste0(metrics$SAMPLE_ID[i], '_',
                                          opt$organism, '_.*_Genes\\.counts$'))
}
names(geneFn) = metrics$SAMPLE_ID

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]

ent_ids = mapIds(org.Rn.eg.db, geneMap$Geneid, "ENTREZID", "ENSEMBL", multiVals='CharacterList')

#  Here we have a potential problem
length(which(sapply(ent_ids, length) > 1)) # 935
length(which(sapply(temp, length) == 2)) # 757

#  For simplicity, we'll check those Entrez IDs which come in a pair for a
#  given Ensembl ID.
ambig_ids = ent_ids[sapply(ent_ids, length) == 2]
first_ids = sapply(ambig_ids, function(x) x[1])
second_ids = sapply(ambig_ids, function(x) x[2])

a = mapIds(org.Rn.eg.db, first_ids, 'GO', 'ENTREZID', multiVals='CharacterList')
b = mapIds(org.Rn.eg.db, second_ids, 'GO', 'ENTREZID', multiVals='CharacterList')

#  Compute a metric to measure similarity in GO results:
#  0 is completely disjoint results; 1 is equivalent results
similarity = rep(0, length(a))
for (i in 1:length(a)) {
    temp = 2 * length(intersect(a[[i]], b[[i]]))
    similarity[i] = temp / (length(setdiff(a[[i]], b[[i]])) + temp)
}

#  This indicates there is a significant difference in GO results for different
#  Entrez IDs matching to one Ensembl ID
mean(similarity) # about 0.49
