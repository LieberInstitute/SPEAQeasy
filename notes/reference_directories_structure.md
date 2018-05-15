### Reference files ###

The pipeline uses many reference files during a run. Due to size limitations in git repositories, not every reference file can be versionated.

The basic Annotation and Genotyping directories are cloned with this repository. On the first run (test run or real run) for a particular species (i. e. hg38, h19, mm10 or rn6), the missing annotation files are built in a species-dependant manner.

Full Annotation and Genotyping directories should looks like this (files buildable by this pipeline are marked):

````
Annotation/
├── chrom_sizes
│   ├── get.chrom.sizes.R
│   ├── hg19.chrom.sizes.gencode
│   ├── hg38.chrom.sizes.gencode
│   ├── mm10.chrom.sizes.gencode
│   └── rn6.chrom.sizes.ensembl
├── ensembl
│   ├── ensemblv75_exon_names_hg19_deduplicated.rda
│   ├── ensemblv75_exon_names_hg19_duplicates_in.rda
│   └── Rnor_6.0
│       ├── fa
│       │   └── Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa (### BUILDABLE BY PIPELINE ###)
│       └── index
│           ├── hisat2_Rnor6.0toplevel.1.ht2 (### BUILDABLE BY PIPELINE ###)
│           ├── hisat2_Rnor6.0toplevel.2.ht2 (### BUILDABLE BY PIPELINE ###)
│           ├── hisat2_Rnor6.0toplevel.3.ht2 (### BUILDABLE BY PIPELINE ###)
│           ├── hisat2_Rnor6.0toplevel.4.ht2 (### BUILDABLE BY PIPELINE ###)
│           ├── hisat2_Rnor6.0toplevel.5.ht2 (### BUILDABLE BY PIPELINE ###)
│           ├── hisat2_Rnor6.0toplevel.6.ht2 (### BUILDABLE BY PIPELINE ###)
│           ├── hisat2_Rnor6.0toplevel.7.ht2 (### BUILDABLE BY PIPELINE ###)
│           └── hisat2_Rnor6.0toplevel.8.ht2 (### BUILDABLE BY PIPELINE ###)
├── ERCC
│   ├── ERCC92.fa
│   └── ERCC92.idx
├── ercc_actual_conc.txt
├── GENCODE
│   ├── exon_names_hg19_hg38_deduplicated.rda
│   ├── exon_names_hg19_hg38.rda
│   ├── GRCh37_hg19
│   │   ├── assembly
│   │   │   ├── fa
│   │   │   │   └── GRCh37.primary_assembly.genome.fa (### BUILDABLE BY PIPELINE ###)
│   │   │   └── index
│   │   │       ├── hisat2_GRCh37primary.1.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh37primary.2.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh37primary.3.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh37primary.4.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh37primary.5.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh37primary.6.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh37primary.7.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       └── hisat2_GRCh37primary.8.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   └── transcripts
│   │       ├── fa
│   │       │   └── gencode.v25lift37.transcripts.fa (### BUILDABLE BY PIPELINE ###)
│   │       └── salmon
│   │           └── salmon_0.8.2_index_gencode.v25lift37.transcripts
│   │               ├── duplicate_clusters.tsv (### BUILDABLE BY PIPELINE ###)
│   │               ├── hash.bin (### BUILDABLE BY PIPELINE ###)
│   │               ├── header.json (### BUILDABLE BY PIPELINE ###)
│   │               ├── indexing.log (### BUILDABLE BY PIPELINE ###)
│   │               ├── quasi_index.log (### BUILDABLE BY PIPELINE ###)
│   │               ├── refInfo.json (### BUILDABLE BY PIPELINE ###)
│   │               ├── rsd.bin (### BUILDABLE BY PIPELINE ###)
│   │               ├── sa.bin (### BUILDABLE BY PIPELINE ###)
│   │               ├── txpInfo.bin (### BUILDABLE BY PIPELINE ###)
│   │               └── versionInfo.json (### BUILDABLE BY PIPELINE ###)
│   ├── GRCh38_hg38
│   │   ├── assembly
│   │   │   ├── fa 
│   │   │   │   └── GRCh38.primary_assembly.genome.fa (### BUILDABLE BY PIPELINE ###)
│   │   │   └── index
│   │   │       ├── hisat2_GRCh38primary.1.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh38primary.2.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh38primary.3.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh38primary.4.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh38primary.5.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh38primary.6.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       ├── hisat2_GRCh38primary.7.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   │       └── hisat2_GRCh38primary.8.ht2 (### BUILDABLE BY PIPELINE ###)
│   │   └── transcripts
│   │       ├── fa
│   │       │   └── gencode.v25.transcripts.fa (### BUILDABLE BY PIPELINE ###)
│   │       └── salmon
│   │           └── salmon_0.8.2_index_gencode.v25.transcripts
│   │               ├── duplicate_clusters.tsv (### BUILDABLE BY PIPELINE ###)
│   │               ├── hash.bin (### BUILDABLE BY PIPELINE ###)
│   │               ├── header.json (### BUILDABLE BY PIPELINE ###)
│   │               ├── indexing.log (### BUILDABLE BY PIPELINE ###)
│   │               ├── quasi_index.log (### BUILDABLE BY PIPELINE ###)
│   │               ├── refInfo.json (### BUILDABLE BY PIPELINE ###)
│   │               ├── rsd.bin (### BUILDABLE BY PIPELINE ###)
│   │               ├── sa.bin (### BUILDABLE BY PIPELINE ###)
│   │               ├── txpInfo.bin (### BUILDABLE BY PIPELINE ###)
│   │               └── versionInfo.json (### BUILDABLE BY PIPELINE ###)
│   └── GRCm38_mm10
│       ├── assembly
│       │   ├── fa
│       │   │   └── GRCm38.primary_assembly.genome.fa (### BUILDABLE BY PIPELINE ###)
│       │   └── index
│       │       ├── GRCm38_mmhisat2_GRCm38primary.1.ht2 (### BUILDABLE BY PIPELINE ###)
│       │       ├── GRCm38_mmhisat2_GRCm38primary.2.ht2 (### BUILDABLE BY PIPELINE ###)
│       │       ├── GRCm38_mmhisat2_GRCm38primary.3.ht2 (### BUILDABLE BY PIPELINE ###)
│       │       ├── GRCm38_mmhisat2_GRCm38primary.4.ht2 (### BUILDABLE BY PIPELINE ###)
│       │       ├── GRCm38_mmhisat2_GRCm38primary.5.ht2 (### BUILDABLE BY PIPELINE ###)
│       │       ├── GRCm38_mmhisat2_GRCm38primary.6.ht2 (### BUILDABLE BY PIPELINE ###)
│       │       ├── GRCm38_mmhisat2_GRCm38primary.7.ht2 (### BUILDABLE BY PIPELINE ###)
│       │       └── GRCm38_mmhisat2_GRCm38primary.8.ht2 (### BUILDABLE BY PIPELINE ###)
│       └── transcripts
│           ├── fa
│           │   └── gencode.vM11.transcripts.fa (### BUILDABLE BY PIPELINE ###)
│           └── salmon
│               └── salmon_0.8.2_index_gencode.vM11.transcripts (### BUILDABLE BY PIPELINE ###)
│                   ├── duplicate_clusters.tsv (### BUILDABLE BY PIPELINE ###)
│                   ├── hash.bin (### BUILDABLE BY PIPELINE ###)
│                   ├── header.json (### BUILDABLE BY PIPELINE ###)
│                   ├── indexing.log (### BUILDABLE BY PIPELINE ###)
│                   ├── quasi_index.log (### BUILDABLE BY PIPELINE ###)
│                   ├── refInfo.json (### BUILDABLE BY PIPELINE ###)
│                   ├── rsd.bin (### BUILDABLE BY PIPELINE ###)
│                   ├── sa.bin (### BUILDABLE BY PIPELINE ###)
│                   ├── txpInfo.bin (### BUILDABLE BY PIPELINE ###)
│                   └── versionInfo.json (### BUILDABLE BY PIPELINE ###)
├── gencode_counts
│   ├── R16-099_Gencode.v25.hg19_Exons.counts
│   ├── R16-099_Gencode.v25.hg19_Genes.counts
│   ├── R16-099_Gencode.v25.hg38_Exons.counts
│   └── R16-099_Gencode.v25.hg38_Genes.counts
├── GENCODE_PRI
│   └── geneCHR.Rdata
├── gs
│   ├── create_gs_gencode_v25_hg19.R
│   ├── create_gs_gencode_v25.R
│   ├── gs_gencode_v25_hg19.Rdata
│   └── gs_gencode_v25_hg38.Rdata
├── junction_txdb
│   ├── feature_to_Tx_ensembl_v75.rda
│   ├── feature_to_Tx_ensembl_v85.rda
│   ├── feature_to_Tx_hg19_gencode_v25lift37.rda
│   ├── feature_to_Tx_hg38_gencode_v25.rda
│   ├── gene_exon_maps_ensembl_v75.rda
│   ├── junction_annotation_hg19_ensembl_v75.rda
│   ├── junction_annotation_hg19_gencode_v25lift37.rda
│   ├── junction_annotation_hg19_refseq_grch37.rda
│   ├── junction_annotation_hg38_ensembl_v85.rda
│   ├── junction_annotation_hg38_gencode_v25.rda
│   ├── junction_annotation_hg38_refseq_grch38.rda
│   ├── junction_annotation_mm10_ensembl_v86.rda
│   ├── junction_annotation_mm10_gencode_vM11.rda
│   └── junction_annotation_rn6_ensembl_v86.rda
└── RSeQC
    ├── hg19
    │   ├── bed
    │   │   └── hg19.bed
    │   └── gtf
    │       └── gencode.v25lift37.annotation.gtf (### BUILDABLE BY PIPELINE ###)
    ├── hg38
    │   ├── bed
    │   │   └── hg38.bed
    │   └── gtf
    │       └── gencode.v25.annotation.gtf (### BUILDABLE BY PIPELINE ###)
    ├── mm10
    │   ├── bed
    │   │   └── mm10.bed
    │   └── gtf
    │       └── gencode.vM11.annotation.gtf (### BUILDABLE BY PIPELINE ###)
    └── rn6
        ├── bed
        │   └── rn6.bed
        └── gtf
            └── Rattus_norvegicus.Rnor_6.0.86.gtf (### BUILDABLE BY PIPELINE ###)
````


````
Genotyping/
├── common_missense_SNVs_hg19.bed
├── common_missense_SNVs_hg38.bed
├── hg19ToHg38.over.chain
├── liftover_variants.R
├── pileVar.pl
└── pull_brain_genotypes.sh
````
