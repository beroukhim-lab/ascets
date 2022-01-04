### ASCETS: Arm-level Somatic Copy-number Events in Targeted Sequencing
### v1.0
### ascets_resources.R

### Author: Liam F. Spurr, liam.spurr@uchospitals.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: October 15, 2020

### License: GNU GPL2, Copyright (C) 2020 Dana-Farber Cancer Institute
### Dependencies: R >= 3.4.0, Libraries: tidyverse, data.table
### See README for guide on how to run this package

#######################################

message("ASCETS: Arm-level Somatic Copy-number Events in Targeted Sequencing")
message("Author: Liam F. Spurr, liam.spurr@uchospitals.edu")
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

###  REQUIRED LIBRARIES
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(data.table)))

### FUNCTIONS

# internal function to manage command line arguments
handle_command_args <- function(args) {
  # make sure all flags are paired (there is an even number) and the required ones have been supplied
  if(length(args) %% 2 != 0 | (!(all(c("-i", "-c", "-o") %in% args)))) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # load in the CNV seg file
  cna <<- as.data.frame(fread(arg_df$value[arg_df$flag == "-i"])) %>% distinct()
  
  # load the remaining required files and parameters
  cytoband <<- as.data.frame(fread(arg_df$value[arg_df$flag == "-c"]))
  min_boc <<- ifelse(length(arg_df$value[arg_df$flag == "-t"]) > 0, as.numeric(arg_df$value[arg_df$flag == "-t"]), 0.5)
  name <<- arg_df$value[arg_df$flag == "-o"]
  
  # check if optional noise file has been supplied
  threshold <<- ifelse(length(arg_df$value[arg_df$flag == "-t"]) > 0, arg_df$value[arg_df$flag == "-t"], 0.2)
  
  # check if the user has specified whether to keep noisy segments (by default noisy segments are excluded)
  keep_noisy <<- ifelse(length(arg_df$value[arg_df$flag == "-k"]) > 0, as.logical(arg_df$value[arg_df$flag == "-k"]), F)
  
  if(length(arg_df$value[arg_df$flag == "-e"]) > 0) noise <<- as.data.frame(fread(arg_df$value[arg_df$flag == "-e"]))
  else noise <<- data.frame()
  
  
  # set amount of arm that must be altered to make an arm-level call
  alteration_threshold <<- ifelse(length(arg_df$value[arg_df$flag == "-a"]) > 0, as.numeric(arg_df$value[arg_df$flag == "-a"]), 0.7)
}

# internal function to make the final alteration call for a chromosome arm
make_arm_call <- function(df, thresh) {
  df <- as.data.frame(df)

  # ensure data frame is not empty
  if(nrow(df) < 1) {return()}
  
  # check to make sure arm passes provided BOC threshold
  if(!df$cov_pass[1]) {return(data.frame(sample = as.character(df$sample[1]), ARM = as.character(df$arm[1]), CALL = "LOWCOV"))}
  
  call <- ""
  
  # determine which alteration type is most prevalent in the arm
  max_alt <- max(df$alt_frac, na.rm = T)
  max_alt_type <- df[df$alt_frac == max_alt,]$alt
  
  # assign call based on alteration threshold
  if(max_alt_type == "NEUTRAL") {
    if(max_alt >= thresh) call <- max_alt_type
    else call <- "NC"
  } else {
    if (max_alt >= thresh) {  # check if the alteration fraction exceeds the threshold for a call
      call <- max_alt_type
    } else {  # otherwise, cannot make a call
      call <- "NC"
    }
  }
  
  # return the result
  data.frame(sample = as.character(df$sample[1]), ARM = as.character(df$arm[1]), CALL = call)
}

# calculates the total area outside the BOC on a chromosome arm
get_missing_cov <- function(df) {
  df <- df %>% arrange(segment_start) # make sure SEG file is properly sorted
  # set counter variable
  tot = 0
  
  # check BOC on the ends of the arm
  if((df[1,3] - df[1,9]) > 0) {tot = tot + (df[1,3] - df[1,9])} # calculate the length from the end of the arm to the first measured base
  if((df[1,10] - df[nrow(df),4]) > 0) {tot = tot + (df[1,10] - df[nrow(df),4])} # calculate the length from the last measured base to the end of the arm
  
  # iterate through middle segments to identify areas outside BOC
  if(nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      diff = df[i,3] - df[i-1,4] # check BOC between consecutive segments
      if (diff > 1) {
        tot = tot + diff # add to the total BOC
      }
    }
  }
  
  # return input data frame with BOC
  df %>% mutate(missing_cov = tot)
}

# general internal function to perform the specified function by chromsome arm for a sample
by_arm <- function(df, f, ...) {
  df_a <- split(df, df$arm) # split the input data frame by arm
  df_out <- lapply(df_a, f, ...) # apply the specified function to each arm
  do.call("rbind", df_out) # provide the result as a data frame
}

# internal wrapper function to get all corrected BOC values for every segment in each sample in a data frame
correct_for_boc <- function(df) {
  df_sam <- split(df, df$sample) # split the input file by sample
  df_sam <- lapply(df_sam, by_arm, get_missing_cov) # correct BOC by arm in each sample
  df_out <- do.call("rbind", df_sam) # provide the result as a data frame
  df_out %>% mutate(cyto_len_corr = cyto_len - missing_cov) # compute the corrected arm length
}

# internal wrapper function to make all aSCNA calls for every sample in a data frame
make_all_calls <- function(df, thresh) {
  df_sam <- split(df, df$sample) # split the input file by sample
  df_sam <- lapply(df_sam, by_arm, make_arm_call, thresh) # make calls for each arm in each sample
  calls <- do.call("rbind", df_sam) # provide the result as a data frame
  calls <- calls %>% mutate(CALL = as.character(CALL)) %>% replace_na(list(CALL = "NC"))
  calls %>% group_by(sample) %>% spread(ARM, CALL, fill = "LOWCOV") %>% ungroup() # turn the output into a matrix
}

# internal function to compute the noise in a given segment
compute_noise <- function(df) {
  if(nrow(df) < 2) return(NA) # requires at least two data points to proceed
  
  # split data into odd and even indexed elements
  seq1 <- df$log2ratio[seq(1, nrow(df), by = 2)]
  seq2 <- df$log2ratio[seq(2, nrow(df), by = 2)]
  
  # compute the absolute difference between the means
  abs(mean(seq1[!is.infinite(seq1)]) - mean(seq2[!is.infinite(seq2)]))
}

# internal function to compute the noise for every segment in a sample
noise_by_segment <- function(df.seg, df.log) {
  # calculate the noise for every segment in the data frame
  n <- lapply(1:nrow(df.seg),
              function(x) compute_noise(df.log %>%
                                          filter(chrom == df.seg$chrom[x],
                                                 start >= df.seg$segment_start[x],
                                                 end <= df.seg$segment_end[x])))
  
  cbind(df.seg, noise = unlist(n)) # return the result
}

# internal function to compute the noise for every sample in a file
noise_by_sample <- function(df.seg, df.log) {
  # calculate the noise for every sample in the data frame
  m <- lapply(unique(df.seg$sample), function(x) noise_by_segment(df.seg %>% filter(sample == !!x),
                                                              df.log %>% filter(sample == !!x)))
  do.call("bind_rows", m) # return the result
}

# internal function to model noise in the segments to determine a noise threshold
determine_noise_threshold <- function(cna, noise, keep_noisy = F) {
  names(noise) <- c("sample", "interval", "gene_exon", "log2ratio")
  missing_from_noise <- length(unique((noise %>% filter(!(sample %in% cna$sample)))$sample))
  if(missing_from_noise > 0) message(paste0(missing_from_noise, " samples in seg file missing from log2 copy ratio file."))
  
  # split the interval column so it can be aligned with the segment file
  noise <- noise %>% separate(interval, into = c("chrom", "start", "end")) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
  
  # calculate the noise in each segment in each sample
  noise <- noise_by_sample(cna %>% filter(sample %in% noise$sample), noise %>% filter(sample %in% cna$sample)) %>% mutate(noise = as.numeric(as.character(noise)))
  
  # determine the threshold at which we are confident that we can make an accurate call (one standard deviation above the mean)
  n <- mean(noise$noise[!is.infinite(noise$noise)], na.rm = T) + sd(noise$noise[!is.infinite(noise$noise)], na.rm = T)
  
  # collect segments where the noise exceeds the threshold
  min_markers <- 5
  noisy <- noise %>% filter(noise > n & num_mark >= min_markers) %>% drop_na()
  
  # calculate the percentage of segments in the file that exceed the thresold
  exceed_frac <- nrow(noisy) * 100 / nrow(cna)
  
  # remove the noisy segments from the segmentation file if specified
  if(!keep_noisy) {
    noisy <- noisy %>% mutate(seg_id = paste0(sample, "_", chrom, ":", segment_start, ":", segment_end))
    cna <- cna %>% mutate(seg_id = paste0(sample, "_", chrom, ":", segment_start, ":", segment_end)) %>%
      filter(!(seg_id %in% noisy$seg_id)) %>% dplyr::select(-seg_id)
  }
  
  # set the amplification and deletion thresholds
  amp_thresh <- n
  del_thresh <- -1 * n
  
  list(cna = cna,
       noise = noise,
       amp_thresh = amp_thresh,
       del_thresh = del_thresh,
       exceed_frac = exceed_frac)
}

# internal function to annotate segments with supplied cytobands
annotate_cytobands <- function(cna, cytoband) {
  left_join(cna, cytoband, by = "chrom") %>% ungroup() %>% filter(segment_start < end & (segment_start > start | segment_end > start)) %>%
    mutate(start = as.numeric(start), end = as.numeric(end),
           segment_start = as.numeric(segment_start), segment_end = as.numeric(segment_end),
           cyto_len = end - start) %>%
    distinct()
}

# internal function to crop segments to bounds of supplied cytobands
crop_segments <- function(cna) {
  cna %>%
    mutate(alt_len = ifelse(end < segment_end & segment_start > start, (end - segment_start), # these correspond to the 4 ways two segments could overlap
                            ifelse(end > segment_end & segment_start > start, (segment_end - segment_start),
                                   ifelse(end < segment_end & segment_start < start, (end - start),
                                          ifelse(end > segment_end & segment_start < start, (segment_end - start), NA)))))
}

# internal function to compute arm weighted average LCRs
compute_weighted_averages <- function(cna) {
  suppressMessages(cna %>% mutate(log2ratio.w = alt_len * log2ratio) %>%
                     group_by(sample, arm) %>%
                     summarize(weight.ave.segmean = round(sum(log2ratio.w) / sum(alt_len), 4)) %>%
                     ungroup() %>%
                     spread(arm, weight.ave.segmean))
}

# internal function to compute arm alteration fractions
compute_alt_fractions <- function(cna, amp_thresh, del_thresh, min_boc) {
  cna %>%
    mutate(alt = ifelse(log2ratio <= del_thresh, "DEL", ifelse(log2ratio >= amp_thresh, "AMP", "NEUTRAL")),
           perc_chrom_cov = cyto_len_corr / cyto_len, # calculate how much of the chromosome was covered
           cov_pass = ifelse(perc_chrom_cov > min_boc, T, F)) %>%
    group_by(sample, arm, cyto_len_corr, alt) %>%
    mutate(alt_len_sum = sum(alt_len)) %>% # get how much of each arm is amplified and deleted
    ungroup() %>%
    mutate(alt_frac = alt_len_sum / cyto_len_corr) %>%  # calculate the fraction
    dplyr::select(sample, arm, alt, perc_chrom_cov, cov_pass, alt_frac) %>%
    distinct()
}

calc_aneu_scores <- function(calls) {
  calls %>% gather(arm, call, -sample) %>%
    group_by(sample) %>%
    summarize(aneuploidy_score = length(call[call %in% c("AMP", "DEL")]) / 
                length(call[call != "LOWCOV"]))
}

# main function to run the ASCETS algorithm

# INPUT
# - CNV segmentation file (*cna* argument supplied as a data frame)
#   - The following columns are required in this order (names can vary): sample, chromosome, segment start, segment end, number of markers, log2ratio
#   - See sample data on Github for an example file
# - Chromosome arm genomic coordinates (*cytoband* argument supplied as a data frame)
#   - See sample data on Github for an example file
# - Minimum arm BOC to make a call (range 0.0 - 1.0; *min_boc* argument supplied as a numeric value) 
#   - Optional, defaults to 0.5
# - Output file prefix (*name* argument supplied as a character string)
# - Individual log2 ratio copy ratios (LCRs) that were used to build copy-number segments (*noise* argument supplied as a data frame)
#   - The following columns are required in this order (names can vary): sample, genomic interval, gene/exon (can be blank), log2ratio
#   - See sample data for an example file
# - A user-defined threshold can also be supplied instead (*threshold* argument supplied as a numeric value)
#   - If this is not supplied, a default segment mean threshold of +/- 0.2 will be used
# - Logical value specifying whether to retain noisy segments (*keep_noisy* argument, suppled as a boolean value: T or F)
#   - Defaults to F
# - Arm level-alteration fraction threshold (*alteration_threshold* argument supplied as a numeric value)
#   - Defaults to 0.7

# OUTPUT:
# A list containing:
# - calls: data frame containing aSCNA calls
# - weight_ave: data frame arm weighted average segment means
# - amp_thresh: numeric value representing amplification threshold that was used to make aSCNA calls
# - del_thresh: numeric value representing deletion threshold that was used to make aSCNA calls
# - alteration_thresh: numeric value representing fraction of arm that must be altered to call an aSCNA
# - name: character string representing name for dataset that was supplied to the function
# - cna: data frame containing segmented CNA data that was supplied to the function
# - noise: data frame containing noise data that was supplied to the function

ascets <- function(cna, cytoband, min_boc = 0.5, name, noise = data.frame(), keep_noisy = F, threshold = 0.2, alteration_threshold = 0.7) {
  
  # enforce names and data classes for seg file
  names(cna) <- c("sample", "chrom", "segment_start", "segment_end", "num_mark", "log2ratio")
  cna <- cna %>% mutate_all(as.character) %>% mutate_at(vars("segment_start", "segment_end", "log2ratio"), function(x) suppressWarnings(as.numeric(x))) %>% drop_na()
  
  ### DETERMINE AMPLIFICATION AND DELETION THRESHOLDS
  if(nrow(noise) > 0) {
    cat("Computing data noise...\n")
    noise_results <- determine_noise_threshold(cna, noise, keep_noisy)
    amp_thresh <- noise_results$amp_thresh
    del_thresh <- noise_results$del_thresh
    cna <- noise_results$cna
    noise <- noise_results$noise
    exceed_frac <- noise_results$exceed_frac
    
  } else {
    amp_thresh <- as.numeric(threshold)
    del_thresh <- -1 * as.numeric(threshold)
    exceed_frac <- NA
    
  }
  
  # align segments to cytobands
  cna_output <- annotate_cytobands(cna, cytoband)
  
  # run BOC correction
  cat("Determining breadth of coverage...\n")
  cna_output <- correct_for_boc(cna_output)
  
  # crop segments to fit provided cytobands
  cat("Cropping segments...\n")
  cna_output <- crop_segments(cna_output)
  
  # save the weighted average segment mean for later output
  cat("Computing weighted averages...\n")
  weight_ave <- compute_weighted_averages(cna_output)
  
  # compute arm alteration fractions
  cat("Computing alteration fractions...\n")
  cna_output <- compute_alt_fractions(cna_output, amp_thresh, del_thresh, min_boc)
  
  # make final arm level calls
  cat("Making arm level calls...\n")
  calls <- make_all_calls(cna_output, alteration_threshold)
  
  # calculate aneuploidy scores
  cat("Calculating aneuploidy scores...\n")
  aneu_scores <- calc_aneu_scores(calls)
  
  list(calls = calls, 
       aneu_scores = aneu_scores,
       weight_ave = weight_ave, 
       amp_thresh = amp_thresh, 
       del_thresh = del_thresh, 
       alteration_thresh = alteration_threshold, 
       exceed_frac = exceed_frac,
       name = name,
       cna = cna,
       noise = noise)
}

# write relevant outputs from ASCETS to a file

# INPUT: a list object output by the ascets() function
# OUTPUT: files are saved to the disk as decribed in the README

write_outputs_to_file <- function(ascets, location = "./") {
  
  if(nrow(ascets$noise) > 0) {
    cat("Outputting noise histogram...\n")
    pdf(paste0(location, ascets$name, "_noise_hist.pdf"))
    hist(ascets$noise$noise, breaks = 100, xlab = "noise", main = "Histogram of noise")
    abline(v = ascets$amp_thresh, lwd = 2, col = "red")
    garbage <- dev.off()
  }
  
  cat("Outputting segment mean histogram...\n")
  pdf(paste0(location, ascets$name, "_segmean_hist.pdf"))
  hist(ascets$cna$log2ratio, breaks = 200, main = "Histogram of LCRs")
  abline(v = ascets$amp_thresh, col = "red", lwd = 2)
  abline(v = ascets$del_thresh, col = "blue", lwd = 2)
  garbage <- dev.off()

  cat("Writing final outputs...\n")
  write.table(ascets$calls, paste0(location, ascets$name, "_arm_level_calls.txt"), quote = F, row.names = F, sep = "\t")
  write.table(ascets$weight_ave, paste0(location, ascets$name, "_arm_weighted_average_segmeans.txt"), quote = F, row.names = F, sep = "\t")
  write.table(ascets$aneu_scores, paste0(location, ascets$name, "_aneuploidy_scores.txt"), quote = F, row.names = F, sep = "\t")
  
  f <- file(paste0(location, ascets$name, "_params.txt"))
  writeLines(c(ascets$name,
               paste0("Amplification LCR threshold: ", ascets$amp_thresh),
               paste0("Deletion LCR threshold: ", ascets$del_thresh),
               paste0("Arm-level alteration call threshold: ", ascets$alteration_thresh),
               paste0("Percent of segments exceeding noise threshold: ", ascets$exceed_frac, "%")), f)
  close(f)
  
  cat("Complete!\n")
}
