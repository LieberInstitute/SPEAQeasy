-- -*- lua -*-
-- vim:ft=lua:et:ts=4

install_dir = ""

help([[
This module loads each of the software tools used by the pipeline.
]])

LmodMessage("Loading RNAsp-software module")

-- ##################################################################
-- Software installed in the repository under .../RNAsp/Software
-- ##################################################################

-- shared folder for several software tools, including bcftools, samtools, and htslib
prepend_path("PATH", install_dir .. "/bin")
prepend_path("MANPATH", install_dir .. "/share/man/man1")


-- fastqc
prepend_path("PATH", install_dir .. "/FastQC")

-- hisat2
prepend_path("PATH", install_dir .. "/hisat2-2.1.0")

-- htslib
prepend_path("C_INCLUDE_PATH", install_dir .. "/include")
prepend_path("LIBRARY_PATH", install_dir .. "/lib")
prepend_path("PKG_CONFIG_PATH", install_dir .. "/lib/pkgconfig")

-- kallisto
prepend_path("PATH", install_dir .. "/kallisto_linux-v0.43.0")

-- regtools
prepend_path("PATH", install_dir .. "/regtools-0.3.0/build")

-- rseqc
prepend_path("PATH", install_dir .. "/RSeQC-2.6.4/bin/usr/local/bin/")
prepend_path("PYTHONPATH", install_dir .. "/RSeQC-2.6.4/lib/")

-- salmon
prepend_path("PATH", install_dir .. "/Salmon-0.8.2_linux_x86_64/bin/")

-- subread (featureCounts)
prepend_path("PATH", install_dir .. "/subread-1.5.0-p3-Linux-x86_64/bin/")

-- trimmomatic
prepend_path("PATH", install_dir .. "/Trimmomatic-0.36/")
prepend_path("PATH", install_dir .. "/Trimmomatic-0.36/adapters/")

-- wiggletools