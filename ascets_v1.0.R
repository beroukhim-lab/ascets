### ASCETS: Arm and chromosome-level Somatic Copy-number Event detection in Targeted Sequencing
### v1.0

### Author: Liam Flinn Spurr, liamf_spurr@dfci.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: July 21, 2020

### License: Pending
### Dependencies: R >= 3.4.0, Libraries: tidyverse, data.table
### See README for guide on how to run this package

#######################################

message("ASCETS: Arm and chromosome-level Somatic Copy-number Event detection in Targeted Sequencing")
message("Author: Liam F. Spurr, liamf_spurr@dfci.harvard.edu")
message("Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu")

### LOAD REQUIRED LIBRARIES
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(data.table)))
###

### FUNCTIONS

# manages command line arguments
handle_command_args <- function(args) {
  # make sure all flags are paired (there is an even number) and the required ones have been supplied
  if(length(args) %% 2 != 0 | length(args) < 8) stop("Command line arguments supplied incorrectly!")

  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)

  # load in the CNV seg file
  cna <<- as.data.frame(fread(arg_df$value[arg_df$flag == "-i"])) %>% distinct()

  # enforce names and data classes
  names(cna) <<- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  cna <<- cna %>% mutate_all(as.character) %>% mutate_at(vars("loc.start", "loc.end", "seg.mean"), function(x) suppressWarnings(as.numeric(x))) %>% drop_na()

  # load the remaining required files and parameters
  cytoband <<- as.data.frame(fread(arg_df$value[arg_df$flag == "-c"]))
  min_cov <<- as.numeric(arg_df$value[arg_df$flag == "-m"])
  prefix <<- arg_df$value[arg_df$flag == "-o"]

  # check if optional noise file has been supplied
  noise_supplied <<- ifelse(length(arg_df$value[arg_df$flag == "-e"]) > 0, T, F)
  threshold_supplied <<- ifelse(length(arg_df$value[arg_df$flag == "-t"]) > 0, T, F)

  if(noise_supplied) {
    # if so load the noise file
    noise <<- as.data.frame(fread(arg_df$value[arg_df$flag == "-e"]))
    # enforce names
    names(noise) <<- c("sample", "interval", "gene_exon", "log2ratio")

    # check if the user has specified whether to keep noisy segments (by default noisy segments are excluded)
    keep_noisy <<- ifelse(length(arg_df$value[arg_df$flag == "-k"]) > 0, as.logical(arg_df$value[arg_df$flag == "-k"]), F)
  }

  if(threshold_supplied) {
    threshold <<- arg_df$value[arg_df$flag == "-t"]
  }
}

# makes the final alteration call for an arm
make_arm_call <- function(df, u, l) {
  df <- as.data.frame(df)

  # ensure data frame is not empty
  if(nrow(df) < 1) {return()}
  call <- ""

  # determine which alteration type is most prevalent in the arm
  max_alt <- max(df$alt_frac)
  max_alt_type <- df[df$alt_frac == max_alt,]$alt

  # assign call based on alteration threshold
  if(max_alt_type == "NEUTRAL") {
    if(max_alt >= u) call <- max_alt_type
    else call <- "NC"
  } else {
    if (max_alt >= u) {  # check if the alteration fraction exceeds the threshold for a call
      call <- max_alt_type
    } else if (max_alt <= l) { # if not, see if it qualifies as a neutral arm
      call <- "NEUTRAL"
    } else {  # otherwise, cannot make a call
      call <- "NC"
    }
  }

  # check to make sure arm passes provided coverage threshold
  if(!df$cov_pass[1]) {call <- "LOWCOV"}

  # return the result
  data.frame(ID = as.character(df$ID[1]), ARM = as.character(df$arm[1]), CALL = call)
}

# gets the missing coverage on an arm
get_missing_cov <- function(df) {
  df <- df %>% arrange(loc.start) # make sure SEG file is properly sorted
  # set coverage counter variable
  tot = 0

  # check if any coverage is missing on the ends of the arm
  if((df[1,3] - df[1,9]) > 0) {tot = tot + (df[1,3] - df[1,9])} # calculate the length from the end of the arm to the first measured base
  if((df[1,10] - df[nrow(df),4]) > 0) {tot = tot + (df[1,10] - df[nrow(df),4])} # calculate the length from the last measured base to the end of the arm

  # iterate through middle segments to identify missing coverage areas
  if(nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      diff = df[i,3] - df[i-1,4] # see if there is any missing coverage between consecutive segments
      if (diff > 1) {
        tot = tot + diff # add any missing coverage to the total
      }
    }
  }

  # return input data frame with length of missing coverage
  df %>% mutate(missing_cov = tot)
}

# general function to perform the specified function by arm for a sample
by_arm <- function(df, f, ...) {
  df_a <- split(df, df$arm) # split the input data frame by arm
  df_out <- lapply(df_a, f, ...) # apply the specified function to each arm
  do.call("rbind", df_out) # provide the result as a data frame
}

# wrapper function to get all corrected coverage values for every sample in a data frame
correct_coverage <- function(df) {
  df_sam <- split(df, df$ID) # split the input file by sample
  df_sam <- lapply(df_sam, by_arm, get_missing_cov) # perform the missing coverage by arm in each sample
  df_out <- do.call("rbind", df_sam) # provide the result as a data frame
  df_out %>% mutate(cyto_len_corr = cyto_len - missing_cov) # compute the corrected arm length
}

# wrapper function to make all arm calls for every sample in a data frame
make_all_calls <- function(df, u, l) {
  df_sam <- split(df, df$ID) # split the input file by sample
  df_sam <- lapply(df_sam, by_arm, make_arm_call, u, l) # make calls for each arm in each sample
  do.call("rbind", df_sam) # provide the result as a data frame
}

# computes the noise in a given segment
compute_noise <- function(df) {
  if(nrow(df) < 2) return(NA) # requires at least two data points to proceed

  # split data into odd and even indexed elements
  seq1 <- df$log2ratio[seq(1, nrow(df), by = 2)]
  seq2 <- df$log2ratio[seq(2, nrow(df), by = 2)]

  # compute the absolute difference between the means
  abs(mean(seq1[!is.infinite(seq1)]) - mean(seq2[!is.infinite(seq2)]))
}

# computes the noise for every segment in a sample
noise_by_segment <- function(df.seg, df.log) {
  # calculate the noise for every segment in the data frame
  n <- lapply(1:nrow(df.seg),
              function(x) compute_noise(df.log %>%
                                          filter(chrom == df.seg$chrom[x],
                                                 start >= df.seg$loc.start[x],
                                                 end <= df.seg$loc.end[x])))

  cbind(df.seg, noise = unlist(n)) # return the result
}

# computes the noise for every sample in a file
noise_by_sample <- function(df.seg, df.log) {
  # calculate the noise for every sample in the data frame
  m <- lapply(unique(df.seg$ID), function(x) noise_by_segment(df.seg %>% filter(ID == !!x),
                                                              df.log %>% filter(sample == !!x)))
  do.call("bind_rows", m) # return the result
}

###


### LOAD INPUT FILES
cat("Parsing data...\n")
handle_command_args(commandArgs(trailingOnly = TRUE))
###

### DETERMINE AMPLIFICATION AND DELETION THRESHOLDS
if(noise_supplied) {
  #if(nrow(noise %>% filter(sample %in% cna$ID)) != 0) stop("Error: samples in noise file do not match seg file")
  cat("Computing data noise...\n")

  missing_from_noise <- length(unique((noise %>% filter(!(sample %in% cna$ID)))$sample))
  if(missing_from_noise > 0) message(paste0(missing_from_noise, " samples in seg file missing from log2 copy ratio file."))

  # split the interval column so it can be aligned with the segment file
  noise <- noise %>% separate(interval, into = c("chrom", "start", "end")) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))

  # calculate the noise in each segment in each sample
  noise <- noise_by_sample(cna %>% filter(ID %in% noise$sample), noise %>% filter(sample %in% cna$ID)) %>% mutate(noise = as.numeric(as.character(noise)))

  # determine the threshold at which we are confident that we can make an accurate call (one standard deviation above the mean)
  n <- mean(noise$noise[!is.infinite(noise$noise)], na.rm = T) + sd(noise$noise[!is.infinite(noise$noise)], na.rm = T)

  # output a histogram of the noise for manual QC
  cat("Outputting noise histogram...\n")
  pdf(paste0(prefix, "_noise_hist.pdf"))
  hist(noise$noise, breaks = 100, xlab = "noise", main = "Histogram of noise")
  abline(v = n, lwd = 2, col = "red")
  garbage <- dev.off()

  # collect segments where the noise exceeds the threshold
  min_markers <- 5
  noisy <- noise %>% filter(noise > n & num.mark >= min_markers) %>% drop_na()

  # calculate the percentage of segments in the file that exceed the thresold
  exceed_frac <- nrow(noisy) * 100 / nrow(cna)

  # remove the noisy segments from the segmentation file if specified
  # these will now be treated as areas of missing coverage
  if(!keep_noisy) {
    noisy <- noisy %>% mutate(seg_id = paste0(ID, "_", chrom, ":", loc.start, ":", loc.end))
    cna <- cna %>% mutate(seg_id = paste0(ID, "_", chrom, ":", loc.start, ":", loc.end)) %>%
      filter(!(seg_id %in% noisy$seg_id)) %>% select(-seg_id)
  }

  # set the amplification and deletion thresholds
  amp_thresh <- n
  del_thresh <- -1 * n

} else if(threshold_supplied) {
  amp_thresh <- as.numeric(threshold)
  del_thresh <- -1 * as.numeric(threshold)
  exceed_frac <- NA

} else {
  # if no noise files are provided, use previously defined thresholds
  amp_thresh <- 0.2
  del_thresh <- -0.2
  exceed_frac <- NA
}

# output histogram distribution of segment mean values for manual QC
cat("Outputting segment mean histogram...\n")
pdf(paste0(prefix, "_segmean_hist.pdf"))
hist(cna$seg.mean, breaks = 200)
abline(v = amp_thresh, col = "red", lwd = 2)
abline(v = del_thresh, col = "blue", lwd = 2)
garbage <- dev.off()
###


### COMPUTE MISSING COVERAGE
# align segments to cytobands
cna_cyto <- left_join(cna, cytoband, by = "chrom") %>% ungroup() %>% filter(loc.start < end & (loc.start > start | loc.end > start)) %>%
  mutate(start = as.numeric(start), end = as.numeric(end),
         loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end),
         cyto_len = end - start) %>%
  distinct()

# run coverage correction
cat("Determining regions of missing coverage...\n")
cna_cyto <- correct_coverage(cna_cyto)

# determine lengths of segments that meet the amplification/deletion thresholds on each arm
cna_cyto <- cna_cyto %>%
  mutate(alt_len = ifelse(end < loc.end & loc.start > start, (end - loc.start), # these correspond to the 4 ways two segments could overlap
                          ifelse(end > loc.end & loc.start > start, (loc.end - loc.start),
                                 ifelse(end < loc.end & loc.start < start, (end - start),
                                        ifelse(end > loc.end & loc.start < start, (loc.end - start), NA)))))

# save the weighted average segment mean for later output
weight_ave <- suppressMessages(cna_cyto %>% mutate(seg.mean.w = alt_len * seg.mean) %>%
  group_by(ID, arm) %>%
  summarize(weight.ave.segmean = round(sum(seg.mean.w) / sum(alt_len), 4)) %>%
  ungroup() %>%
  spread(arm, weight.ave.segmean))

cna_cyto <- cna_cyto %>%
  mutate(alt = ifelse(seg.mean <= del_thresh, "DEL", ifelse(seg.mean >= amp_thresh, "AMP", "NEUTRAL")), # only retain segments that meet the thresholds
         perc_chrom_cov = cyto_len_corr / cyto_len, # calculate how much of the chromosome was covered
         cov_pass = ifelse(perc_chrom_cov > min_cov, T, F)) %>%
  group_by(ID, arm, cyto_len_corr, alt) %>%
  mutate(alt_len_sum = sum(alt_len)) %>% # get how much of each arm is amplified and deleted
  ungroup() %>%
  mutate(alt_frac = alt_len_sum / cyto_len_corr) %>%  # calculate the fraction
  select(ID, arm, alt, perc_chrom_cov, cov_pass, alt_frac) %>%
  distinct()

# define alteration (upper) and neutral (lower) thresholds
upper_thresh <- 0.7
lower_thresh <- 0.3

# output a histogram of alteration fractions for manual QC
cat("Outputting alteration fraction histogram...\n")
pdf(paste0(prefix, "_altfrac_hist.pdf"))
hist(cna_cyto$alt_frac[cna_cyto$alt != "NEUTRAL"],
     xlab = "Fraction of arm altered (stratified by amp/del)",
     main = "Histogram of cohort-wide arm alterations")
abline(v = lower_thresh, col = "#B24C63", lwd = 2)
abline(v = upper_thresh, col = "#23CE6B", lwd = 2)
garbage <- dev.off()
###

### MAKE ARM LEVEL CALLS
# make final arm level calls
cat("Making arm level calls...\n")
calls <- make_all_calls(cna_cyto, upper_thresh, lower_thresh)
calls <- calls %>% mutate(CALL = as.character(CALL)) %>% replace_na(list(CALL = "NC"))
calls_out <- calls %>% group_by(ID) %>% spread(ARM, CALL, fill = "LOWCOV") %>% ungroup() # turn the output into a matrix

# write output to a file
cat("Writing output...\n")
write.table(calls_out, paste0(prefix, "_arm_level_calls.txt"), quote = F, row.names = F, sep = "\t")
write.table(weight_ave, paste0(prefix, "_arm_weighted_average_segmeans.txt"), quote = F, row.names = F, sep = "\t")

# write algorithm parameters to a file
f <- file(paste0(prefix, "_params.txt"))
writeLines(c(prefix,
             paste0("Amplification segment mean threshold: ", amp_thresh),
             paste0("Deletion segment mean threshold: ", del_thresh),
             paste0("Arm-level alteration call threshold: ", upper_thresh),
             paste0("Neutral call threshold: ", lower_thresh),
             paste0("Percent of segments exceeding noise threshold: ", exceed_frac, "%")), f)
close(f)

cat("Complete!\n")
###
