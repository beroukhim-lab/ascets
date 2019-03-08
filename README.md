# ACCENTS
**A**rm and **C**hromosome **C**opy number **E**vents through **N**oise modeling in **T**argeted **S**equencing

This repository contains the code to run the ACCENTS arm and chromosome-level copy number events caller for targeted sequencing data.

## How to run
### System requirements
This code has been run on Mac and Linux using R version 3.4.0. The user must also install R packages "tidyverse" and 
"data.table" from CRAN prior to running the program.

### Input files

Required:
- CNV segmentation file (.seg)
	- We use .seg files from RobustCNV in our pipeline
- Chromosome arm genomic coordinates
	- We supply an example file in the repository for hg19, but other coordinates can be supplied in the same format
- Minimum arm coverage to make a call (range 0.0 - 1.0)
	- We usually use 0.5
	- Specify 0.0 to allow any coverage
- Output file prefix

Optional (recommended):
- log2 ratios that were used to build segments
	- We use exonic calls from RobustCNV
	- If this is not supplied, a previously-defined segment mean threshold of +/- 0.2 will be used
- Logical value specifying whether to remove noisy segments (T or F)
	- Noisy segments will be treated as no coverage areas

### Run
```bash
Rscript accents_v1.0.R \
	-i [CNV segment file] \ 
	-c [genomic chromosome arm coordinates] \ 
	-m [minimum arm coverage] \
	-n [exonic log2 ratios used to compute noise] [OPTIONAL] \
	-k [logical specifying whether to remove noisy segments] [OPTIONAL] \
	-o [file output prefix]
```
