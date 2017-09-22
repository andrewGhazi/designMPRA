library(tidyverse)
library(pwr)
library(parallel)

# For a given data file
computeTransfectionStatistics = function(fileNum){
  # Get the file's depth normalization factor
  depthNum = fileDepths %>% filter(src == file_names[fileNum]) %>% .$fileDepth 
  
  # Read in the counts
  rnaCounts = read_tsv(paste0(dir, file_names[fileNum]), 
                       col_names = c('allele', 'barcode', 'count')) %>% 
    mutate(depthAdjCount = 1e6*count/depthNum) #normalize them for depth
  
  # Then compute statistics for the file
  alleleStatistics = depthAdjDNAmeanCount %>% 
    filter(bcMean > .06) %>% # this is where we introduce the cutoff 
    left_join(rnaCounts, by = c('allele', 'barcode')) %>% #join onto DNA counts
    mutate(activity = log(depthAdjCount/bcMean)) %>% #compute activity
    group_by(allele) %>% #for each allele
    summarise(alleleMean = mean(activity, na.rm = TRUE), #compute statistics
              alleleSD = sd(activity, na.rm = TRUE),
              numBarcodes = sum(!is.na(count)),
              lillieforsP = ifelse(numBarcodes > 4, 
                                   lillie.test(na.omit(activity))$p.value, 
                                   NA),
              n = n()) %>% 
    filter(numBarcodes > 1) %>% # take only alleles which had > 1 barcode 
    mutate(file = file_names[fileNum]) # and add on a file identifier
  # The filter removes alleles that only had zero or one barcodes show up in
  # the RNA, for which a standard deviation is not meaningful
  
  return(alleleStatistics)
}

get_nonnormal_allele = function(file_allele_list){
  dir = '/mnt/bigData2/andrew/MPRA/Tewhey/indivTags/'
  
  file_name = file_allele_list[['file_name']]
  file_name_str = file_name
  nonnormal_allele = file_allele_list[['allele']]
  
  depthNum = fileDepths %>% filter(src == paste0('Geuv_90K_', file_name, '.tag.ct.indiv.expanded')) %>% .$fileDepth 
  
  rnaCounts = read_tsv(paste0(dir, 'Geuv_90K_', file_name, '.tag.ct.indiv.expanded'), 
                       col_names = c('allele', 'barcode', 'count')) %>% 
    mutate(depthAdjCount = 1e6*count/depthNum) #normalize them for depth
  
  # return a data_frame with the activity measurements for one allele
  depthAdjDNAmeanCount %>% 
    filter(bcMean > .06) %>% # this is where we introduce the cutoff 
    left_join(rnaCounts, by = c('allele', 'barcode')) %>% #join onto DNA counts
    na.omit %>% 
    filter(allele == nonnormal_allele) %>% 
    mutate(activity = log(depthAdjCount/bcMean),
           file_name = file_name_str) # and add on a file identifier
  # The filter removes alleles that only had zero or one barcodes show up in
  # the RNA, for which a standard deviation is not meaningful
}

load('~/designMPRA/outputs/tewheyDepthAdjDNAcounts.RData')
load('~/designMPRA/outputs/tewheyAllTransfectionsStats.RData')

dir = '/mnt/bigData2/andrew/MPRA/Tewhey/indivTags/'

file_names = list.files(dir, 
                        pattern = '.expanded$')

getFileDepth = function(file){
  # function to get the total number of reads in a sample. 
  # Used to normalize counts by sample depth.
  
  read_tsv(paste0(dir, file),
           col_names = c('allele', 'barcode', 'count')) %>% 
    .$count %>% 
    sum
}

#apply the above function to the sample files
fileDepths = data_frame(src = file_names,
                        fileDepth = mclapply(src, 
                                             getFileDepth,
                                             mc.cores = 6) %>% unlist)

# worst = list(allele = 'rs28519456_RC_altB', file_name = 'HepG2.r5') %>%
#   get_nonnormal_allele()

good = list(allele = 'rs12942988_altB', file_name = 'HepG2.r5') %>% 
  get_nonnormal_allele()

good_mean = mean(good$activity)

# used to calculate effect size
pooled_SD = sqrt(((nrow(good) - 1)*sd(good$activity)**2 + (nrow(good) - 1)*.926**2) / (nrow(good)*2 - 2))

simulate_power = function(transcription_shift, sig_level = .05, n_sim = 10000){
  
  #simulate a t-test many times
  monte_carlo_p_values = sapply(1:n_sim, function(x){
    
    # randomly sample the most non-normal allele
    non_normal_samples = sample(good$activity, size = nrow(good), replace = TRUE) 
    
    # randomly sample a "typical" allele of the same size
    normal_samples = rnorm(nrow(good), 
                           mean = good_mean - transcription_shift, 
                           sd = .926) # using the average SD from Ulirsch and Tewhey together
    
    # compute the p-value of a t-test
    t.test(non_normal_samples, normal_samples)$p.value
  })
  
  # return the fraction that pass the input significance level
  sum(monte_carlo_p_values < sig_level) / n_sim
}

sim_power = data_frame(transcription_shift = seq(-1.1, 1.1, length.out = 50),
                       theoretical_power = pwr.t.test(n = nrow(good),
                                                      d = abs(transcription_shift) / pooled_SD,
                                                      type = 'two.sample',
                                                      alternative = 'two.sided',
                                                      sig.level = 1e-5)$power,
                       monte_carlo_power = mclapply(transcription_shift, 
                                                    simulate_power,
                                                    mc.cores = 8,
                                                    sig_level = 1e-5, n_sim = 3e6) %>% unlist)

save(good, sim_power, file = '~/designMPRA/outputs/sim_power_rs12942988_altB.RData')