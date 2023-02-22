#   Pull FASTQ files used to generate the SPEAQeasy test data. Note that for
#   human, paired-end, forward data, the SPEAQeasy files
#   "dm3_file[12]_[12].fastq" came from
#   https://github.com/nellore/rail/tree/master/ex files
#   "dm3_example_[12]_(left|right).fastq", respectively, and the files
#   "sample_0[12]_[12].fastq.gz" come from
#   https://github.com/LieberInstitute/RNAseq-pipeline/tree/master/ex_simulated/paired_end_stranded.

#  Human paired reverse
wget https://sra-pub-src-1.s3.amazonaws.com/SRR9600560/HFNYCBBXX_102879-001-005_CTTAGGAC_L001_R1.fastq.1
wget https://sra-pub-src-1.s3.amazonaws.com/SRR9600560/HFNYCBBXX_102879-001-005_CTTAGGAC_L001_R2.fastq.1

wget https://sra-pub-src-1.s3.amazonaws.com/SRR9600558/HFNYCBBXX_102879-001-006_ATCTGACC_L001_R1.fastq.1
wget https://sra-pub-src-1.s3.amazonaws.com/SRR9600558/HFNYCBBXX_102879-001-006_ATCTGACC_L001_R2.fastq.1

#  Human single reverse

wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12424559/LungOrganoid-10uM-Imatinib-3_R1.fastq.gz.1
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12424558/LungOrganoid-10uM-Imatinib-2_R1.fastq.gz.1

#  Human single forward

wget https://sra-pub-src-2.s3.amazonaws.com/ERR2812363/human_testis_rna_3.fastq.gz.1
wget https://sra-pub-src-2.s3.amazonaws.com/ERR2812362/human_testis_rna_2.fastq.gz.1

#  Mouse paired reverse

wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12405003/D14_3_1.fq.gz.1
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12405003/D14_3_2.fq.gz.1

wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12405004/D14_2_1.fq.gz.1
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12405004/D14_2_2.fq.gz.1

#  Mouse single reverse

wget https://sra-download.ncbi.nlm.nih.gov/traces/sra68/SRZ/012456/SRR12456219/sham_rep2_R1_fastq.gz

wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR337/ERR3374068/278-4.fastq.gz


#  Mouse single forward (we have this, but in case we want a known source)

wget https://sra-pub-src-2.s3.amazonaws.com/SRR10700194/Clozapine_Colon_4.fastq.gz.1

wget https://sra-pub-src-2.s3.amazonaws.com/ERR2812389/mouse_liver_ribo_2.fastq.gz.1

wget https://sra-pub-src-2.s3.amazonaws.com/SRR10700255/Dox2mgmL_Lung_3.fastq.gz.1

#  Rat paired unstranded (not needed)

wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12663720/siPRMT5_3.R1.fastq.gz.1
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12663720/siPRMT5_3.R2.fastq.gz.1

wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12321291/FXST5_1.fq.gz.1
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR12321291/FXST5_2.fq.gz.1

#  Rat paired reverse

wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR11680484/CTW0_20171012_AGATGT_S34_L003_R1_001.fastq.gz
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR11680484/CTW0_20171012_AGATGT_S34_L003_R2_001.fastq.gz

wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR11680485/CTW1_20171012_AGCACC_S35_L003_R1_001.fastq.gz
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR11680485/CTW1_20171012_AGCACC_S35_L003_R2_001.fastq.gz

#  Rat single forward (not needed)

wget https://sra-pub-src-2.s3.amazonaws.com/SRR9825112/AB3010.fastq.gz.1
#wget https://sra-pub-src-2.s3.amazonaws.com/SRR9825113/AB3011.fastq.gz.1

wget https://sra-download.ncbi.nlm.nih.gov/traces/sra46/SRZ/011821/SRR11821109/Day_2_Hepatocyte_4.1.HT5CMBGXC_4.2.fq.gz

#  Rat single reverse

wget https://sra-download.ncbi.nlm.nih.gov/traces/sra64/SRZ/011486/SRR11486935/REF2c_S10_L001_R1_001.fastq.gz
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra61/SRZ/011486/SRR11486936/REF2c_S10_L002_R1_001.fastq.gz
