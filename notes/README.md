# Development Notes for the LIEBER INSTITUTE JAFFE LAB RNA - SEQ ANALYSIS PIPELINE #

## Completed functionalities pending release ##

	+ Complete functionality for **paired-end** type of data for **mouse (mm10)**.

		+ Variant Calling 
		+ Expressed Regions detection
		+ Full Coverage Rdata generation
		+ Transcript Counts Rdata generation

	+ Functional End processes for **paired-end** type of data for **human (hg19, hg38)**.

		+ Variant Calling
		+ Expressed Regions detection
		+ Full Coverage Rdata generation

	+ Functional End processes for **single-end** AND **paired-end** type of data for **rat (rn6)**.

		+ Variant Calling
		+ Expressed Regions detection

## Known Issues ##

	+ **Human** (hg19, hg38) **paired-end**

		+ Transcript Counts Rdata generation:
			+ The R script `create_count_objects-human.R` has some troubles reading FastQC data files.
			+ When the mouse version `create_count_objects-mouse.R` was adapted to read human files, the process was successful.
			+ Thus, some minor bugs must persist in the human version script. A possible solution is to transplant the PE code block in the mouse script, if Lieber team approves.

	+ **rat** (hg19, hg38) **paired-end** and **single-end**

		+ Full Coverage Rdata generation
			+ Bugs persist. Fixing requires testing and debugging in the Lieber Environment.
		+ Transcript Counts Rdata generation
			+ Bugs persist. Fixing requires testing and debugging in the Lieber Environment.

**Process validation status**.

  + For runs in any mode (System, Docker, and SGE)
 
 Green: validated process; Yellow: process with a pinpointed bug; Red: Process with issues not pinpointed hitherto; Gray: Process not used by that run.
 
![Validations](https://github.com/LieberInstitute/RNAsp/blob/feature/Code_revision/notes/Process_Validation_table.png)