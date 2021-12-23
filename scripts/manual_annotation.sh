#  This script exists for users who do not have access to an internet
#  connection on the machine(s) where the pipeline will be executed. Run this
#  script inside the SPEAQeasy repository directory on a machine with internet
#  access to prepare the pipeline for runs without internet access. The
#  necessary annotation files are pulled/built. See configuration below.

#############################################################
#  User configuration
#############################################################

reference="hg38"
config="slurm.config"
annotation_dir="$PWD/Annotation" # This should match the "--annotation" flag in your main script
work_dir="$PWD/manual_annotation_work" # Location of temporary working directory

#  This should be "docker", "singularity", or "local"
installation_method="singularity"


#############################################################
#  Version-dependent variables
#############################################################

r_container="libddocker/bioc_kallisto:3.13"
r_version="4.1.0"

#############################################################
#  Determine links and filenames
#############################################################

#  Verify validity of installation method
if [[ ! ("$installation_method" == "docker" || "$installation_method" == "singularity" || "$installation_method" == "local") ]]; then
    echo "The 'installation_method' variable must either be 'docker', 'singularity', or 'local'."
    exit 1
fi

#  Grab the annotation versions specified at the above config
gencode_version_human=$(cat conf/$config | grep "gencode_version_human" | cut -d "=" -f 2 | tr -d ' |"')
gencode_version_mouse=$(cat conf/$config | grep "gencode_version_mouse" | cut -d "=" -f 2 | tr -d ' |"')
ensembl_version_rat=$(cat conf/$config | grep "ensembl_version_rat" | cut -d "=" -f 2 | tr -d ' |"')
anno_build=$(cat conf/$config | grep "anno_build" | cut -d "=" -f 2 | cut -d " " -f 2 | tr -d '"')

origDir=$PWD
mkdir -p $work_dir

if [ $reference == "hg38" ]; then
    anno_version=$gencode_version_human
    anno_suffix="${reference}_gencode_v${anno_version}_$anno_build"
    num_chrs=25
    
    # Reference assembly fasta, gtf, and transcript fasta
    fa_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${anno_version}/GRCh38.primary_assembly.genome.fa.gz"
    if [ $anno_build == "primary" ]; then
        gtf_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${anno_version}/gencode.v${anno_version}.primary_assembly.annotation.gtf.gz"
    else
        gtf_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${anno_version}/gencode.v${anno_version}.annotation.gtf.gz"
    fi
    tx_fa_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${anno_version}/gencode.v${anno_version}.transcripts.fa.gz"
    
elif [ $reference == "hg19" ]; then
    echo "Warning: use of 'primary' annotation is not supported for hg19, as GENCODE does not provide a primary .gtf file. Continuing with annotation build 'main'."
    anno_build="main"
        
    anno_version=$gencode_version_human
    anno_suffix="${reference}_gencode_v${anno_version}lift37_${anno_build}"
    num_chrs=25
    
    # Reference assembly fasta, gtf, and transcript fasta
    fa_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${anno_version}/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
    gtf_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${anno_version}/GRCh37_mapping/gencode.v${anno_version}lift37.annotation.gtf.gz"
    tx_fa_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${anno_version}/GRCh37_mapping/gencode.v${anno_version}lift37.transcripts.fa.gz"
    
elif [ $reference == "mm10" ]; then
    anno_version=$gencode_version_mouse
    anno_suffix="${reference}_gencode_${anno_version}_$anno_build"
    num_chrs=22
    
    # Reference assembly fasta, gtf, and transcript fasta
    fa_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${anno_version}/GRCm38.primary_assembly.genome.fa.gz"
    if [ $anno_build == "primary" ]; then
        gtf_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${anno_version}/gencode.v${anno_version}.primary_assembly.annotation.gtf.gz"
    else
        gtf_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${anno_version}/gencode.v${anno_version}.annotation.gtf.gz"
    fi
    tx_fa_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${anno_version}/gencode.v${anno_version}.transcripts.fa.gz"
    
else  # rat
    anno_version=$ensembl_version_rat
    anno_suffix="${reference}_ensembl_${anno_version}_$anno_build"
    num_chrs=23
    
    # Reference assembly fasta, gtf, and transcript fasta
    fa_link="ftp://ftp.ensembl.org/pub/release-${anno_version}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
    if [ $anno_build == "primary" ]; then
        gtf_link="ftp://ftp.ensembl.org/pub/release-${anno_version}/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.${anno_version}.gtf.gz"
    else
        gtf_link="ftp://ftp.ensembl.org/pub/release-${anno_version}/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.${anno_version}.chr.gtf.gz"
    fi
    tx_fa_link="ftp://ftp.ensembl.org/pub/release-${anno_version}/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"
fi

#############################################################
#  Download the primary assembly fasta and subset to
#  build the corresponding "main" fasta
#############################################################

destDir="${annotation_dir}/reference/${reference}/assembly/fa"
mkdir -p $destDir
cd $destDir

#  Name of the primary assembly fasta after being downloaded and unzipped
baseFile=$(basename $(echo $fa_link | rev | cut -d "." -f 2- | rev))
            
#  Name the pipeline will use for the primary and main assembly fastas, respectively
primaryFile="assembly_$(echo $anno_suffix | rev | cut -d "_" -f 2- | rev)_primary.fa"
mainFile="assembly_$(echo $anno_suffix | rev | cut -d "_" -f 2- | rev)_main.fa"

#  Name of fasta to use for this pipeline execution instance
out_fasta="assembly_${anno_suffix}.fa"

#  ----------------------------
#  Pull primary assembly fasta
#  ----------------------------
echo "Pulling and unzipping assembly fasta..."
curl -O $fa_link
gunzip ${baseFile}.gz
mv $baseFile $primaryFile # rename for consistency with pipeline naming conventions

#  ----------------------------
#  Build main assembly fasta
#  ----------------------------
echo "Subsetting primary fasta to build main fasta..."

#  Find the line of the header for the first extra contig (to not
#  include in the "main" annotation fasta
first_bad_line=$(grep -n ">" $primaryFile | cut -d : -f 1 | paste -s | cut -f $(($num_chrs + 1)))
            
#  Make a new file out of all the lines up and not including that
sed -n "1,$(($first_bad_line - 1))p;${first_bad_line}q" $primaryFile > $mainFile

cp $out_fasta $work_dir/


#############################################################
#  Download the gtf
#############################################################

destDir="${annotation_dir}/RSeQC/${reference}/gtf"
mkdir -p $destDir
cd $destDir

#  Names of gtf file when downloaded + unzipped and after renamed, respectively
baseFile=$(basename $(echo $gtf_link | rev | cut -d "." -f 2- | rev))
out_gtf="transcripts_${anno_suffix}.gtf"

#  Pull, unzip, and rename gtf
echo "Pulling and unzipping the gtf file..."
curl -O $gtf_link
gunzip ${baseFile}.gz
mv $baseFile $out_gtf

#  Copy the gtf to the work directory
cp $out_gtf $work_dir/

#############################################################
#  Download the transcripts fasta
#############################################################

destDir="${annotation_dir}/reference/${reference}/transcripts/fa"
mkdir -p $destDir
cd $destDir

#  Name of transcripts fasta when downloaded
baseFile=$(basename $(echo $tx_fa_link | rev | cut -d "." -f 2- | rev))

echo "Pulling and unzipping the transcripts fasta..."
curl -O $tx_fa_link
gunzip ${baseFile}.gz

#############################################################

cd $work_dir

if [[ "$installation_method" == "docker" ]; then
    echo "User wants to use docker; building annotation objects inside the R container..."
    cp $origDir/scripts/build_annotation_objects.R $work_dir/
    
    docker run \
        -v $work_dir/:/$work_dir/ \
        $r_container Rscript /$work_dir/build_annotation_objects.R \
            -r $reference \
            -s $anno_suffix
    
    echo "Copying objects out of the container to their destinations..."
    docker cp ${r_container}:/chrom_sizes_${anno_suffix} ${annotation_dir}/junction_txdb/
    docker cp ${r_container}:/junction_annotation_${anno_suffix}.rda ${annotation_dir}/junction_txdb/
    if [ ! $reference == "rn6" ]; then
        docker cp ${r_container}:/feature_to_Tx_${anno_suffix}.rda ${annotation_dir}/junction_txdb/
    fi
elif [[ "$installation_method" == "singularity" ]; then
    echo "Building annotation objects inside the singularity R container..."
    cp $origDir/scripts/build_annotation_objects.R $work_dir/
    
    #  Grab the filename for the R docker image as it exists in the cache
    #  (it is placed in the cache during installation w/ 'singularity' mode
    image_name=$(echo $r_container | sed 's/[:\/]/-/g').img
    
    singularity exec \
        --pwd $PWD \
        -B $work_dir/:/$work_dir/ \
        $origDir/dockerfiles/singularity_cache/$image_name Rscript /$work_dir/build_annotation_objects.R \
            -r $reference \
            -s $anno_suffix
    
    #  Move the generated annotation files to their final destinations
    mv chrom_sizes_${anno_suffix} ${annotation_dir}/junction_txdb/
    mv junction_annotation_${anno_suffix}.rda ${annotation_dir}/junction_txdb/
    if [ ! $reference == "rn6" ]; then
        mv feature_to_Tx_${anno_suffix}.rda ${annotation_dir}/junction_txdb/
    fi
else
    echo "User is using locally installed software; building annotation objects using local R install..."
    $origDir/Software/R-${r_version}/bin/Rscript $origDir/scripts/build_annotation_objects.R \
        -r $reference \
        -s $anno_suffix
    
    echo "Moving annotation objects to their destinations..."
    mv chrom_sizes_${anno_suffix} ${annotation_dir}/junction_txdb/
    mv junction_annotation_${anno_suffix}.rda ${annotation_dir}/junction_txdb/
    if [ ! $reference == "rn6" ]; then
        mv feature_to_Tx_${anno_suffix}.rda ${annotation_dir}/junction_txdb/
    fi
fi

echo "Done configuring annotation. You may now run the pipeline with or without an internet connection."
