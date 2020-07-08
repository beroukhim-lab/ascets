# ACCENTS
**A**rm and **C**hromosome **C**opy-number **E**vent detection in **N**eoplastic **T**argeted **S**equencing

This repository contains the code to run the ACCENTS arm and chromosome-level copy number events caller for targeted sequencing data. ACCENTS produces arm-level copy-number variant calls and arm-level weighted average log2 segment means from segmented copy number data.

## How to run
### System requirements
This code has been run on Mac and Linux using R version 3.4.0. The user must also install R packages "tidyverse" and
"data.table" from CRAN prior to running the program.

### Input files

Required:
- CNV segmentation file (.seg)
	- We use .seg files from RobustCNV in our pipeline
- Chromosome arm genomic coordinates
	- We supply an example file in the repository for hg19 (original data from [bioMart](http://grch37.ensembl.org/biomart/martview/69a5479f5796c22ca786f81386e2d5e4)), but other coordinates can be supplied in the same format such as cytoband coordinates (also provided for hg19)
- Minimum arm coverage to make a call (range 0.0 - 1.0)
	- We usually use 0.5
	- Specify 0.0 to allow any level of coverage
- Output file prefix

Optional (recommended):
- Exonic log2 ratio measurements that were used to build copy-number segments
	- We use exonic calls from RobustCNV
	- If this is not supplied, a previously-defined segment mean threshold of +/- 0.2 will be used
- Logical value specifying whether to remove noisy segments (T or F)
	- Noisy segments will be treated as no coverage areas

### Run (command-line)
```bash
Rscript accents_v1.2.R \
	-i [CNV segment file] \
	-c [genomic chromosome arm coordinates (or other coordinates)] \
	-m [minimum arm coverage] \
	-e [exonic log2 ratios used to compute noise] [OPTIONAL] \
	-k [logical specifying whether to keep noisy segments] [OPTIONAL] \
	-t [manual log2 ratio threshold to determine amplifications and deletions] [OPTIONAL, will overwrite noise estimate] \
	-o [file output prefix]
```

### Output files
- Histogram of the alteration fraction of all arms in the input cohort (altfrac_hist.pdf)
- Arm-level calls for each arm in each sample in the input (arm_level_calls.txt)
- Weighted average segment mean values for each arm in each sample in the input (weighted_average_segmeans.txt)
- Histogram of modeled noise in the segments in the input cohort (noise_hist.pdf)
- Parameters used to make arm-level calls (params.txt)
- Histogram of the segment means in the input cohort (segmean_hist.pdf)
