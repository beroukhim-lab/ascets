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

### LOAD REQUIRED LIBRARIES
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(data.table)))
###

getCurrentFileLocation <-  function() {
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) ==0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

### LOAD INPUT FILES
cat("Loading required functions...\n")
source(paste0(getCurrentFileLocation(), "/ascets_resources.R"))

cat("Parsing data...\n")
handle_command_args(commandArgs(trailingOnly = TRUE))
ascets_output <- ascets(cna, cytoband, min_boc, name, noise, keep_noisy, threshold, alteration_threshold)
write_outputs_to_file(ascets_output)
