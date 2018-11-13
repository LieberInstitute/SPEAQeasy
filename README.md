# WinterFlow Sample

This feature branch is a simplified version of the original RNAsp pipeline.

It features one process, the hisat2 index building, developed under the modularity framework.

To test this proposal, run the following: `nextflow run new_ref_builder.nf --reference "hg38"`

The command will download the human genome fasta from GENCODE (release 29), and use it to build a hisat2 index.

# WinterFlow features

* Self-contained code modules using the mk command tool to establish and follow dependency rules between files (inputs - intermediates - outputs).

* Easier Develop, Testing and debugging for complex pipelines such as Lieber's RNAsp.

* Better, more detailed documentation at code level, but also as sub-readme files in markdown language.

## Contact
iaguilaror@gmail.com
