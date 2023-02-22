#   Subset the SPEAQeasy test data to a small number of reads
N=1000000

#  Human paired reverse
head -n $N FASTQ_raw/HFNYCBBXX_102879-001-005_CTTAGGAC_L001_R1.fastq.1 | gzip -c > FASTQ_subsetted/SRR9600560_R1.fastq.gz
head -n $N FASTQ_raw/HFNYCBBXX_102879-001-005_CTTAGGAC_L001_R2.fastq.1 | gzip -c > FASTQ_subsetted/SRR9600560_R2.fastq.gz

head -n $N FASTQ_raw/HFNYCBBXX_102879-001-006_ATCTGACC_L001_R1.fastq.1 | gzip -c > FASTQ_subsetted/SRR9600558_R1.fastq.gz
head -n $N FASTQ_raw/HFNYCBBXX_102879-001-006_ATCTGACC_L001_R2.fastq.1 | gzip -c > FASTQ_subsetted/SRR9600558_R2.fastq.gz

#  Human single reverse

gunzip -c FASTQ_raw/LungOrganoid-10uM-Imatinib-3_R1.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/human/single/lung3.fastq.gz
gunzip -c FASTQ_raw/LungOrganoid-10uM-Imatinib-2_R1.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/human/single/lung2.fastq.gz

#  Human single forward

gunzip -c FASTQ_raw/human_testis_rna_3.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/human/single/test3.fastq.gz
gunzip -c FASTQ_raw/human_testis_rna_2.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/human/single/test2.fastq.gz

#  Mouse paired reverse

gunzip -c FASTQ_raw/D14_3_1.fq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/mouse/paired/D14_3_1.fq.gz
gunzip -c FASTQ_raw/D14_3_2.fq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/mouse/paired/D14_3_2.fq.gz

gunzip -c FASTQ_raw/D14_2_1.fq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/mouse/paired/D14_2_1.fq.gz
gunzip -c FASTQ_raw/D14_2_2.fq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/mouse/paired/D14_2_2.fq.gz

#  Mouse single reverse

gunzip -c FASTQ_raw/sham_rep2_R1_fastq.gz | head -n $N | gzip -c > FASTQ_subsetted/mouse/single/sham_rep2.fastq.gz

gunzip -c FASTQ_raw/278-4.fastq.gz | head -n $N | gzip -c > FASTQ_subsetted/mouse/single/sample4.fastq.gz

#  Mouse single forward (we have this, but in case we want a known source)

gunzip -c FASTQ_raw/Clozapine_Colon_4.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/mouse/single/Clozapine_Colon_4.fastq.gz

gunzip -c FASTQ_raw/mouse_liver_ribo_2.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/mouse/single/liver_ribo.fastq.gz

gunzip -c FASTQ_raw/Dox2mgmL_Lung_3.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/mouse/single/lung3.fastq.gz

#  Rat paired unstranded (not needed)

gunzip -c FASTQ_raw/siPRMT5_3.R1.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/rat/paired/si_1.fastq.gz
gunzip -c FASTQ_raw/siPRMT5_3.R2.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/rat/paired/si_2.fastq.gz

gunzip -c FASTQ_raw/FXST5_1.fq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/rat/paired/fx_1.fastq.gz
gunzip -c FASTQ_raw/FXST5_2.fq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/rat/paired/fx_2.fastq.gz

#  Rat paired reverse

gunzip -c FASTQ_raw/CTW0_20171012_AGATGT_S34_L003_R1_001.fastq.gz | head -n $N | gzip -c > FASTQ_subsetted/rat/paired/ct_1.fastq.gz
gunzip -c FASTQ_raw/CTW0_20171012_AGATGT_S34_L003_R2_001.fastq.gz | head -n $N | gzip -c > FASTQ_subsetted/rat/paired/ct_2.fastq.gz

gunzip -c FASTQ_raw/CTW1_20171012_AGCACC_S35_L003_R1_001.fastq.gz | head -n $N | gzip -c > FASTQ_subsetted/rat/paired/w1_1.fastq.gz
gunzip -c FASTQ_raw/CTW1_20171012_AGCACC_S35_L003_R2_001.fastq.gz | head -n $N | gzip -c > FASTQ_subsetted/rat/paired/w1_2.fastq.gz

#  Rat single forward (not needed)

gunzip -c FASTQ_raw/AB3010.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/rat/single/ab10.fastq.gz
#gunzip -c FASTQ_raw/AB3011.fastq.gz.1 | head -n $N | gzip -c > FASTQ_subsetted/rat/single/ab11.fastq.gz

gunzip -c FASTQ_raw/Day_2_Hepatocyte_4.1.HT5CMBGXC_4.2.fq.gz | head -n $N | gzip -c > FASTQ_subsetted/rat/single/h4.fastq.gz

#  Rat single reverse

gunzip -c FASTQ_raw/REF2c_S10_L001_R1_001.fastq.gz | head -n $N | gzip -c > FASTQ_subsetted/rat/single/r1.fastq.gz
gunzip -c FASTQ_raw/REF2c_S10_L002_R1_001.fastq.gz | head -n $N | gzip -c > FASTQ_subsetted/rat/single/r2.fastq.gz
