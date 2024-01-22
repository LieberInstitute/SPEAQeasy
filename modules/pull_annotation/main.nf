// This module contains processes related to pulling reference genomes,
// transcriptomes, and annotation files. All such processes are only run once
// for a particular set of annotation settings, and cached locally for future
// SPEAQeasy runs

// Determine the URLS to the reference genome, rerference transcriptome, and
// GTF of annotation for the particular species and 'anno_build' being used
if (params.custom_anno == "") {
    if (params.reference == "hg38") {
        params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh38.primary_assembly.genome.fa.gz"
        if (params.anno_build == "primary") {
            params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/gencode.v${params.anno_version}.primary_assembly.annotation.gtf.gz"
        } else {
            params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/gencode.v${params.anno_version}.annotation.gtf.gz"
        }
        params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/gencode.v${params.anno_version}.transcripts.fa.gz"
    } else if (params.reference == "hg19") {
        params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
        params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/gencode.v${params.anno_version}lift37.annotation.gtf.gz"
        params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/gencode.v${params.anno_version}lift37.transcripts.fa.gz"
    } else if (params.reference == "mm10") {
        params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/GRCm38.primary_assembly.genome.fa.gz"
        if (params.anno_build == "primary") {
            params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/gencode.v${params.anno_version}.primary_assembly.annotation.gtf.gz"
        } else {
            params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/gencode.v${params.anno_version}.annotation.gtf.gz"
        }
        params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/gencode.v${params.anno_version}.transcripts.fa.gz"
    } else { // rat
        // At Ensembl, the rat genome switches from "Rnor_6.0" to "mRatBN7.2" at
        // and after release 105
        if (Integer.parseInt(params.anno_version) >= 105) {
            rat_genome_name = "mRatBN7.2" 
        } else {
            rat_genome_name = "Rnor_6.0" 
        }

        params.fa_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.${rat_genome_name}.dna.toplevel.fa.gz"
        if (params.anno_build == "primary") {
            params.gtf_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/gtf/rattus_norvegicus/Rattus_norvegicus.${rat_genome_name}.${params.anno_version}.gtf.gz"
        } else {
            params.gtf_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/gtf/rattus_norvegicus/Rattus_norvegicus.${rat_genome_name}.${params.anno_version}.chr.gtf.gz"
        }
        params.tx_fa_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.${rat_genome_name}.cdna.all.fa.gz"
    }
}

/*
 * Pull assembly fasta (from GENCODE for human/ mouse, or ensembl for rat).
 * Build a "main" fasta of only canonical sequences from the pulled "primary"
 * fasta containing all sequences/contigs.
 */

// If any files are not already downloaded/ prepared: download the primary
// assembly fasta and create a subsetted fasta of "main" (reference) sequences only
process PullAssemblyFasta {

    tag "Downloading Assembly FA File: ${baseName}"
    storeDir "${params.annotation}/reference/${params.reference}/assembly/fa"
        
    output:
        path "${out_fasta}", emit: reference_fasta
        path "*.fa" // to store the primary and main fastas, regardless of which is used

    shell:
        //  Name of the primary assembly fasta after being downloaded and unzipped
        baseName = file("${params.fa_link}").getName() - ".gz"
        
        //  Name the pipeline will use for the primary and main assembly fastas, respectively
        primaryName = "assembly_${params.anno_suffix}.fa".replaceAll("main", "primary")
        mainName = "assembly_${params.anno_suffix}.fa".replaceAll("primary", "main")
        
        //  Name of fasta to use for this pipeline execution instance
        out_fasta = "assembly_${params.anno_suffix}.fa"
        
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Pull and unzip primary assembly fasta
        curl -O "!{params.fa_link}"
        gunzip "!{baseName}.gz"
        mv !{baseName} !{primaryName} # rename for consistency with pipeline naming conventions
        
        #######################################################################
        #  Create the "main" fasta of canonical seqs only
        #######################################################################
        
        #  Determine how many chromosomes/seqs to keep
        if [ !{params.reference_type} == "human" ]; then
            num_chrs=25
        elif [ !{params.reference} == "mm10" ]; then
            num_chrs=22
        else
            num_chrs=23
        fi
        
        #  Find the line of the header for the first extra contig (to not
        #  include in the "main" annotation fasta
        first_bad_line=$(grep -n ">" !{primaryName} | cut -d : -f 1 | paste -s | cut -f $(($num_chrs + 1)))
        
        #  Make a new file out of all the lines up and not including that
        sed -n "1,$(($first_bad_line - 1))p;${first_bad_line}q" !{primaryName} > !{mainName}
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}


/*
 * Download reference .gtf (from GENCODE for human/ mouse, ensembl for rat)
 */
// Uses "storeDir" to download gtf only when it doesn't exist, and output the cached
// file if it does already exist
process PullGtf {

    tag "Downloading GTF File: ${baseName}"
    storeDir "${params.annotation}/RSeQC/${params.reference}/gtf"

    output:
        path "${out_gtf}", emit: reference_gtf

    shell:
        // Names of gtf file when downloaded + unzipped and after renamed, respectively
        baseName = file("${params.gtf_link}").getName() - ".gz"
        out_gtf = "transcripts_${params.anno_suffix}.gtf"
        '''
        #  Pull, unzip, and rename transcript gtf
        curl -O "!{params.gtf_link}"
        gunzip "!{baseName}.gz"
        mv !{baseName} !{out_gtf}
        '''
}

/*
 * Download transcript FASTA
 */
			
// Uses "storeDir" to download files only when they don't exist, and output the cached
// files if they do already exist
process PullTranscriptFasta {
            
    tag "Downloading TX FA File: ${baseName}"
    storeDir "${params.annotation}/reference/${params.reference}/transcripts/fa"
    
    input:
        path subset_script

    output:
        path baseName, emit: transcript_fa

    shell:
        //  For human and mouse, only the "main" transcripts FASTA is
        //  available. These means the output FASTA won't differ between
        //  "main" and "primary" runs for these species. For rat, only a
        //  "primary" FASTA is available, so the output will differ
        if (params.reference == 'rat') {
            baseName = "transcripts_${params.anno_suffix}.fa"
        } else {
            baseName = file("${params.tx_fa_link}").getName() - ".gz"
        }

        '''
        curl -o !{baseName}.gz !{params.tx_fa_link}
        gunzip !{baseName}.gz

        #   For rat, the only transcripts FASTA available includes "primary"
        #   transcripts. Subset if the user selects "main" build
        if [[ (!{params.reference} == 'rat') && (!{params.anno_build} == "main") ]]; then
            !{params.Rscript} !{subset_script}
        fi
        '''
}
