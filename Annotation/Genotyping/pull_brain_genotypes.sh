#### VCF
plink --bfile /dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10 --extract <(cut -f4 /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed) --recode vcf --out LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_GenotypingBarcode
bgzip LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_GenotypingBarcode.vcf
tabix -p vcf LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_GenotypingBarcode.vcf.gz

## A-transpose
plink --bfile /dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10 --extract <(cut -f4 /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed) --recode A-transpose --out LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_GenotypingBarcode
