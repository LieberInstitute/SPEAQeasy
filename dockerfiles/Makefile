#### Repository Info #### 

NS = libddocker

######################### 
######### Images ########
## R Base
R_IMAGE = r3.4.3_base
R_IMAGE_TAG = 1
R_IMAGE_VERSION = 2
## Ubuntu Base
UBUNTU_BASE_IMAGE = ubuntu16.04_base
UBUNTU_BASE_IMAGE_TAG = 1
UBUNTU_BASE_IMAGE_VERSION = 2
## Kallisto
ERCC_IMAGE = kallisto_v0.43.1
ERCC_IMAGE_TAG = 1
ERCC_IMAGE_VERSION = 2
## FastQC
QUALITY_IMAGE = fastqc_v0.11.5
QUALITY_IMAGE_TAG = 1
QUALITY_IMAGE_VERSION = 2
## Trimmomatic
TRIM_IMAGE = trimmomatic-0.36
TRIM_IMAGE_TAG = 1
TRIM_IMAGE_VERSION = 2
## HISAT 
HISAT_IMAGE = hisat2-2.0.4
HISAT_IMAGE_TAG = 1
HISAT_IMAGE_VERSION = 2
## RSeQC
RSEQC_IMAGE = rseqc-v2.6.4
RSEQC_IMAGE_TAG = 1
RSEQC_IMAGE_VERSION = 2
## Samtools
SAMTOOLS_IMAGE = samtools-1.3.1
SAMTOOLS_IMAGE_TAG = 1
SAMTOOLS_IMAGE_VERSION = 2
## Salmon
SALMON_IMAGE = salmon-0.9.1
SALMON_IMAGE_TAG = 1
SALMON_IMAGE_VERSION = 2
## Regtools
REGTOOLS_IMAGE = regtools-0.3.0
REGTOOLS_IMAGE_TAG = 1
REGTOOLS_IMAGE_VERSION = 2
## SubRead 
SUBREAD_IMAGE = subread-1.6.0
SUBREAD_IMAGE_TAG = 1
SUBREAD_IMAGE_VERSION = 2
## Wiggletools
WIGGLETOOLS_IMAGE = wiggletools-1.2
WIGGLETOOLS_IMAGE_TAG = 1
WIGGLETOOLS_IMAGE_VERSION = 2
## BCFTools
BCFTOOLS_IMAGE = bcftools-1.3.1
BCFTOOLS_IMAGE_TAG = 1
BCFTOOLS_IMAGE_VERSION = 2
## VCFTOols
VCFTOOLS_IMAGE = vcftools-v0.1.15
VCFTOOLS_IMAGE_TAG = 1
VCFTOOLS_IMAGE_VERSION = 2


#########################
###### Image Build ######

build:
	# R Image Base
	docker build -f $(R_IMAGE).dockerfile -t $(NS)/$(R_IMAGE) .
	docker tag $(NS)/$(R_IMAGE) $(NS)/$(R_IMAGE):$(R_IMAGE_TAG)_v$(R_IMAGE_VERSION)
	docker tag $(NS)/$(R_IMAGE) $(NS)/$(R_IMAGE):latest
	# Ubuntu Image Base
	docker build -f $(UBUNTU_BASE_IMAGE).dockerfile -t $(NS)/$(UBUNTU_BASE_IMAGE) .
	docker tag $(NS)/$(UBUNTU_BASE_IMAGE) $(NS)/$(UBUNTU_BASE_IMAGE):$(UBUNTU_BASE_IMAGE_TAG)_v$(UBUNTU_BASE_IMAGE_VERSION)
	docker tag $(NS)/$(UBUNTU_BASE_IMAGE) $(NS)/$(UBUNTU_BASE_IMAGE):latest
	# Kallisto
	docker build -f $(ERCC_IMAGE).dockerfile -t $(NS)/$(ERCC_IMAGE) .
	docker tag $(NS)/$(ERCC_IMAGE) $(NS)/$(ERCC_IMAGE):$(ERCC_IMAGE_TAG)_v$(ERCC_IMAGE_VERSION)
	docker tag $(NS)/$(ERCC_IMAGE) $(NS)/$(ERCC_IMAGE):latest
	# FastQC
	docker build -f $(QUALITY_IMAGE).dockerfile -t $(NS)/$(QUALITY_IMAGE) .
	docker tag $(NS)/$(QUALITY_IMAGE) $(NS)/$(QUALITY_IMAGE):$(QUALITY_IMAGE_TAG)_v$(QUALITY_IMAGE_VERSION)
	docker tag $(NS)/$(QUALITY_IMAGE) $(NS)/$(QUALITY_IMAGE):latest
	# Trimmomatic
	docker build -f $(TRIM_IMAGE).dockerfile -t $(NS)/$(TRIM_IMAGE) .
	docker tag $(NS)/$(TRIM_IMAGE) $(NS)/$(TRIM_IMAGE):$(TRIM_IMAGE_TAG)_v$(TRIM_IMAGE_VERSION)
	docker tag $(NS)/$(TRIM_IMAGE) $(NS)/$(TRIM_IMAGE):latest
	# Hisat
	docker build -f $(HISAT_IMAGE).dockerfile -t $(NS)/$(HISAT_IMAGE) .
	docker tag $(NS)/$(HISAT_IMAGE) $(NS)/$(HISAT_IMAGE):$(HISAT_IMAGE_TAG)_v$(HISAT_IMAGE_VERSION)
	docker tag $(NS)/$(HISAT_IMAGE) $(NS)/$(HISAT_IMAGE):latest
	# RSeQC
	docker build -f $(RSEQC_IMAGE).dockerfile -t $(NS)/$(RSEQC_IMAGE) .
	docker tag $(NS)/$(RSEQC_IMAGE) $(NS)/$(RSEQC_IMAGE):$(RSEQC_IMAGE_TAG)_v$(RSEQC_IMAGE_VERSION)
	docker tag $(NS)/$(RSEQC_IMAGE) $(NS)/$(RSEQC_IMAGE):latest
	# Samtools
	docker build -f $(SAMTOOLS_IMAGE).dockerfile -t $(NS)/$(SAMTOOLS_IMAGE) .
	docker tag $(NS)/$(SAMTOOLS_IMAGE) $(NS)/$(SAMTOOLS_IMAGE):$(SAMTOOLS_IMAGE_TAG)_v$(HISAT_IMAGE_VERSION)
	docker tag $(NS)/$(SAMTOOLS_IMAGE) $(NS)/$(SAMTOOLS_IMAGE):latest
	# Salmon
	docker build -f $(SALMON_IMAGE).dockerfile -t $(NS)/$(SALMON_IMAGE) .
	docker tag $(NS)/$(SALMON_IMAGE) $(NS)/$(SALMON_IMAGE):$(SALMON_IMAGE_TAG)_v$(SALMON_IMAGE_VERSION)
	docker tag $(NS)/$(SALMON_IMAGE) $(NS)/$(SALMON_IMAGE):latest
	# Regtools
	docker build -f $(REGTOOLS_IMAGE).dockerfile -t $(NS)/$(REGTOOLS_IMAGE) .
	docker tag $(NS)/$(REGTOOLS_IMAGE) $(NS)/$(REGTOOLS_IMAGE):$(REGTOOLS_IMAGE_TAG)_v$(REGTOOLS_IMAGE_VERSION)
	docker tag $(NS)/$(REGTOOLS_IMAGE) $(NS)/$(REGTOOLS_IMAGE):latest
	# Subread
	docker build -f $(SUBREAD_IMAGE).dockerfile -t $(NS)/$(SUBREAD_IMAGE) .
	docker tag $(NS)/$(SUBREAD_IMAGE) $(NS)/$(SUBREAD_IMAGE):$(SUBREAD_IMAGE_TAG)_v$(SUBREAD_IMAGE_VERSION)
	docker tag $(NS)/$(SUBREAD_IMAGE) $(NS)/$(SUBREAD_IMAGE):latest
	# Wiggletools
	docker build -f $(WIGGLETOOLS_IMAGE).dockerfile -t $(NS)/$(WIGGLETOOLS_IMAGE) .
	docker tag $(NS)/$(WIGGLETOOLS_IMAGE) $(NS)/$(WIGGLETOOLS_IMAGE):$(WIGGLETOOLS_IMAGE_TAG)_v$(WIGGLETOOLS_IMAGE_VERSION)
	docker tag $(NS)/$(WIGGLETOOLS_IMAGE) $(NS)/$(WIGGLETOOLS_IMAGE):latest
	# BCFTools
	docker build -f $(BCFTOOLS_IMAGE).dockerfile -t $(NS)/$(BCFTOOLS_IMAGE) .
	docker tag $(NS)/$(BCFTOOLS_IMAGE) $(NS)/$(BCFTOOLS_IMAGE):$(BCFTOOLS_IMAGE_TAG)_v$(BCFTOOLS_IMAGE_VERSION)
	docker tag $(NS)/$(BCFTOOLS_IMAGE) $(NS)/$(BCFTOOLS_IMAGE):latest
	# VCFTools
	docker build -f $(VCFTOOLS_IMAGE).dockerfile -t $(NS)/$(VCFTOOLS_IMAGE) .
	docker tag $(NS)/$(VCFTOOLS_IMAGE) $(NS)/$(VCFTOOLS_IMAGE):$(VCFTOOLS_IMAGE_TAG)_v$(VCFTOOLS_IMAGE_VERSION)
	docker tag $(NS)/$(VCFTOOLS_IMAGE) $(NS)/$(VCFTOOLS_IMAGE):latest

#######################
##### Image Push ######

push:
	# R Image
	docker push $(NS)/$(R_IMAGE):$(R_IMAGE_TAG)_v$(R_IMAGE_VERSION)
	docker push $(NS)/$(R_IMAGE):latest
	# Ubuntu Image
	docker push $(NS)/$(UBUNTU_BASE_IMAGE):$(UBUNTU_BASE_IMAGE_TAG)_v$(UBUNTU_BASE_IMAGE_VERSION)
	docker push $(NS)/$(UBUNTU_BASE_IMAGE):latest
	# Kallisto
	docker push $(NS)/$(ERCC_IMAGE):$(ERCC_IMAGE_TAG)_v$(ERCC_IMAGE_VERSION)
	docker push $(NS)/$(ERCC_IMAGE):latest
	# FastQC
	docker push $(NS)/$(QUALITY_IMAGE):$(QUALITY_IMAGE_TAG)_v$(QUALITY_IMAGE_VERSION)
	docker push $(NS)/$(QUALITY_IMAGE):latest
	# Trimmomatic
	docker push $(NS)/$(TRIM_IMAGE):$(TRIM_IMAGE_TAG)_v$(TRIM_IMAGE_VERSION)
	docker push $(NS)/$(TRIM_IMAGE):latest
	# Hisat
	docker push $(NS)/$(HISAT_IMAGE):$(HISAT_IMAGE_TAG)_v$(HISAT_IMAGE_VERSION)
	docker push $(NS)/$(HISAT_IMAGE):latest
	# RSeQC
	docker push $(NS)/$(RSEQC_IMAGE):$(RSEQC_IMAGE_TAG)_v$(RSEQC_IMAGE_VERSION)
	docker push $(NS)/$(RSEQC_IMAGE):latest
	# Samtools
	docker push $(NS)/$(SAMTOOLS_IMAGE):$(SAMTOOLS_IMAGE_TAG)_v$(SAMTOOLS_IMAGE_VERSION)
	docker push $(NS)/$(SAMTOOLS_IMAGE):latest
	# Salmon
	docker push $(NS)/$(SALMON_IMAGE):$(SALMON_IMAGE_TAG)_v$(SALMON_IMAGE_VERSION)
	docker push $(NS)/$(SALMON_IMAGE):latest
	# Regtools
	docker push $(NS)/$(REGTOOLS_IMAGE):$(REGTOOLS_IMAGE_TAG)_v$(REGTOOLS_IMAGE_VERSION)
	docker push $(NS)/$(REGTOOLS_IMAGE):latest
	# Subread
	docker push $(NS)/$(SUBREAD_IMAGE):$(SUBREAD_IMAGE_TAG)_v$(SUBREAD_IMAGE_VERSION)
	docker push $(NS)/$(SUBREAD_IMAGE):latest
	# Wiggletools
	docker push $(NS)/$(WIGGLETOOLS_IMAGE):$(WIGGLETOOLS_IMAGE_TAG)_v$(WIGGLETOOLS_IMAGE_VERSION)
	docker push $(NS)/$(WIGGLETOOLS_IMAGE):latest
	# BCFTools
	docker push $(NS)/$(BCFTOOLS_IMAGE):$(BCFTOOLS_IMAGE_TAG)_v$(BCFTOOLS_IMAGE_VERSION)
	docker push $(NS)/$(BCFTOOLS_IMAGE):latest
	# VCFTools
	docker push $(NS)/$(VCFTOOLS_IMAGE):$(VCFTOOLS_IMAGE_TAG)_v$(VCFTOOLS_IMAGE_VERSION)
	docker push $(NS)/$(VCFTOOLS_IMAGE):latest

#######################
##### Image Depoy #####

deploy: build push

#######################
##### Image Pull ######

pull:
	# R Image
	docker pull $(NS)/$(R_IMAGE):$(R_IMAGE_TAG)_v$(R_IMAGE_VERSION)
	# Ubuntu Image
	docker pull $(NS)/$(UBUNTU_BASE_IMAGE):$(UBUNTU_BASE_IMAGE_TAG)_v$(UBUNTU_BASE_IMAGE_VERSION)
	# Kallisto
	docker pull $(NS)/$(ERCC_IMAGE):$(ERCC_IMAGE_TAG)_v$(ERCC_IMAGE_VERSION)
	# FastQC
	docker pull $(NS)/$(QUALITY_IMAGE):$(QUALITY_IMAGE_TAG)_v$(QUALITY_IMAGE_VERSION)
	# Trimmomatic
	docker pull $(NS)/$(TRIM_IMAGE):$(TRIM_IMAGE_TAG)_v$(TRIM_IMAGE_VERSION)
	# Hisat
	docker pull $(NS)/$(HISAT_IMAGE):$(HISAT_IMAGE_TAG)_v$(HISAT_IMAGE_VERSION)
	# RSeQC
	docker pull $(NS)/$(RSEQC_IMAGE):$(RSEQC_IMAGE_TAG)_v$(RSEQC_IMAGE_VERSION)
	# Samtools
	docker pull $(NS)/$(SAMTOOLS_IMAGE):$(SAMTOOLS_IMAGE_TAG)_v$(SAMTOOLS_IMAGE_VERSION)
	# Salmon
	docker pull $(NS)/$(SALMON_IMAGE):$(SALMON_IMAGE_TAG)_v$(SALMON_IMAGE_VERSION)
	# Regtools
	docker pull $(NS)/$(REGTOOLS_IMAGE):$(REGTOOLS_IMAGE_TAG)_v$(REGTOOLS_IMAGE_VERSION)
	# Subread
	docker pull $(NS)/$(SUBREAD_IMAGE):$(SUBREAD_IMAGE_TAG)_v$(SUBREAD_IMAGE_VERSION)
	# Wiggletools
	docker pull $(NS)/$(WIGGLETOOLS_IMAGE):$(WIGGLETOOLS_IMAGE_TAG)_v$(WIGGLETOOLS_IMAGE_VERSION)
	# BCFTools
	docker pull $(NS)/$(BCFTOOLS_IMAGE):$(BCFTOOLS_IMAGE_TAG)_v$(BCFTOOLS_IMAGE_VERSION)
	# VCFTools
	docker pull $(NS)/$(VCFTOOLS_IMAGE):$(VCFTOOLS_IMAGE_TAG)_v$(VCFTOOLS_IMAGE_VERSION)
