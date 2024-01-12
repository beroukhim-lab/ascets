# ASCETS
**A**rm-level **S**omatic **C**opy-number **E**vents in **T**argeted **S**equencing

_Copyright (C) 2020 Dana-Farber Cancer Institute_

Author: Liam F. Spurr, liam.spurr@uchospitals.edu

Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu

_Please cite: Spurr LF, Touat M, Taylor AM, et al. Quantification of aneuploidy in targeted sequencing data using ASCETS, Bioinformatics (2020), https://doi.org/10.1093/bioinformatics/btaa980._

This repository contains the code to run the ASCETS arm-level copy number events caller for targeted sequencing data. ASCETS produces arm-level copy-number variant calls and arm-level weighted average log2 segment means from segmented copy number data.

## How to run ASCETS
### System requirements
This code has been tested on Mac and Linux using R version 3.6.1. The user must also install R packages [tidyverse](https://www.tidyverse.org/packages/) and
[data.table](https://github.com/Rdatatable/data.table/wiki/Installation) from CRAN prior to running the program.

### Command line

ASCETS can be run directly from the command line to facilitate easy and effecient generation of aSCNA calls from user data.

#### Input files

Required:
- CNV segmentation file (.seg)
	- We use .seg files from RobustCNV in our pipeline
	- The following columns are required in this order (names can vary): sample, chromosome, segment start, segment end, number of markers, log2ratio
	- See sample data for an example file
- Chromosome arm genomic coordinates
	- We supply an example file in the repository for hg19 (original data from [UCSC](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz)), but other coordinates can be supplied in the same format such as cytoband coordinates (also provided for hg19)
- Minimum arm breadth of coverage (BOC) to make a call (range 0.0 - 1.0)
	- Optional, defaults to 0.5
	- Specify 0.0 to allow any BOC
- Output file prefix

Optional (recommended):
- Arm level-alteration fraction threshold
	- Defaults to 0.7
- Individual log2 ratio copy ratios (LCRs) that were used to build copy-number segments
	- We use exonic calls from RobustCNV
	- The following columns are required in this order (names can vary): sample, genomic interval, gene/exon (can be blank), log2ratio
	- See sample data for an example file
	- If this is not supplied, a default segment mean threshold of +/- 0.2 will be used
	- A user-defined threshold can also be supplied instead
- Logical value specifying whether to retain noisy segments (T or F)
	- Defaults to F

#### Run the script

Note: must be run from the same directory as the *run_ascets.R* and *ascets_resources.R* files
```bash
Rscript run_ascets.R \
	-i [CNV segment file] \
	-c [genomic chromosome arm coordinates (or other coordinates)] \
	-m [minimum arm BOC] [OPTIONAL] \
	-a [arm alteration fraction threshold] [OPTIONAL] \
	-e [LCRs used to compute noise] [OPTIONAL] \
	-k [logical specifying whether to keep noisy segments] [OPTIONAL] \
	-t [manual log2 ratio threshold to determine amplifications and deletions] [OPTIONAL, will overwrite noise estimate] \
	-o [file output prefix]
```

#### Sample command
```bash
Rscript run_ascets.R -i seg_example.seg -c genomic_arm_coordinates_hg19.txt -m 0.5 -e lcr_example.txt -k F -a 0.7 -o sample_output
```

#### Output files
- Arm-level calls for each arm in each sample in the input (arm_level_calls.txt)
- Weighted average segment mean values for each arm in each sample in the input (weighted_average_segmeans.txt)
- Aneuploidy scores (range 0-1, fraction of arms amplified or deleted out of total called arms [call ≠ LOWCOV])
- Histogram of modeled noise in the segments in the input cohort (noise_hist.pdf)
- Parameters used to make arm-level calls (params.txt)
- Histogram of the segment means in the input cohort (segmean_hist.pdf)


### R Studio

ASCETS can also be run from within R Studio through calling the *ascets()* function from another script. All required functions are included in the *ascets_resources.R* file. These functions can be imported and ASCETS can be run using the commands below. Please note that the required data files must be supplied with the column order shown in the sample data. Data may be written to a file using the *write_outputs_to_file()* function.

#### Input files

- CNV segmentation file (*cna* argument supplied as a data frame)
	- The following columns are required in this order (names can vary): sample, chromosome, segment start, segment end, number of markers, log2ratio
	- See sample data for an example file
- Chromosome arm genomic coordinates (*cytoband* argument supplied as a data frame)
	- We supply a example files in the repository for hg19/38 (original data from [bioMart](http://grch37.ensembl.org/biomart/martview/69a5479f5796c22ca786f81386e2d5e4)), but other coordinates can be supplied in the same format such as cytoband coordinates (also provided for hg19)
- Minimum arm breadth of coverage (BOC) to make a call (range 0.0 - 1.0; *min_cov* argument supplied as a numeric value) 
	- Optional, defaults to 0.5
	- Specify 0.0 to allow any BOC
- Output file prefix (*name* argument supplied as a character string)
- Individual log2 ratio copy ratios (LCRs) that were used to build copy-number segments (*noise* argument supplied as a data frame)
	- The following columns are required in this order (names can vary): sample, genomic interval, gene/exon (can be blank), log2ratio
	- See sample data for an example file
	- A user-defined threshold can also be supplied instead (*threshold* argument supplied as a numeric value)
	- If this is not supplied, a default segment mean threshold of +/- 0.2 will be used
- Logical value specifying whether to retain noisy segments (*keep_noisy* argument, suppled as a boolean value: T or F)
	- Defaults to F
- Arm level-alteration fraction threshold (*alteration_threshold* argument supplied as a numeric value)
	- Defaults to 0.7

#### Output files

A list containing:
- *calls*: data frame containing aSCNA calls
- *weight_ave*: data frame containing arm weighted average segment means
- *aneu_scores*: data frame containing aneuploidy scores for each sample
- *amp_thresh*: numeric value representing amplification threshold that was used to make aSCNA calls
- *del_thresh*: numeric value representing deletion threshold that was used to make aSCNA calls
- *alteration_thresh*: numeric value representing fraction of arm that must be altered to call an aSCNA
- *name*: character string representing name for dataset that was supplied to the function
- *cna*: data frame containing segmented CNA data that was supplied to the function
- *noise*: data frame containing noise data that was supplied to the function

#### Run the script

```r
source("ascets_resources.R")
ascets_output <- ascets(cna, cytoband, min_cov, name, noise, keep_noisy, threshold, alteration_threshold)
write_outputs_to_file(ascets_output, location = "output_folder/")
```
