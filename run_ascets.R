### ASCETS: Arm-level Somatic Copy-number Events in Targeted Sequencing
### v1.0
### run_ascets.R

### Author: Liam F, Spurr, liam.spurr@uchospitals.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: October 15, 2020

### License: GNU GPL2, Copyright (C) 2020 Dana-Farber Cancer Institute
### Dependencies: R >= 3.4.0, Libraries: tidyverse, data.table
### See README for guide on how to run this package

#######################################

message("ASCETS: Arm-level Somatic Copy-number Events in Targeted Sequencing")
message("Author: Liam F. Spurr, liamf_spurr@dfci.harvard.edu")
message("Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu\n")

message("Copyright (C) 2020 Dana-Farber Cancer Institute")
message("This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.\n")

### LOAD REQUIRED LIBRARIES
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(data.table)))
###

### LOAD INPUT FILES
cat("Loading required functions...\n")
source("./ascets_resources.R")

cat("Parsing data...\n")
handle_command_args(commandArgs(trailingOnly = TRUE))
ascets_output <- ascets(cna, cytoband, min_cov, name, noise, keep_noisy, threshold, alteration_threshold)
write_outputs_to_file(ascets_output)
