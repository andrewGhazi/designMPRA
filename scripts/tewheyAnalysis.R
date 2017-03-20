# Let's try to see what the Tewhey data is like

library(tidyverse)
library(magrittr)
library(DESeq2)

tew = read_tsv('/mnt/bigData2/andrew/MPRA/Tewhey/GSE75661_79k_collapsed_counts.txt') %>% 
  separate(Oligo, into = c('snp', 'allele'), sep = '(?=[AB]$)') %>% # dat regex #_(?=[A-Za-z]+$)
  gather(block, count, Plasmid_r1:HepG2_r5) %>% 
  mutate(alt = grepl('alt', allele))

getActivities = function(snpDat){
  snpDat %>% group_by(replicate) %>% 
    summarise(Aact = A[grepl('NA', block)] / A[grepl('Plas', block)],
              Bact = B[grepl('NA', block)] / B[grepl('Plas', block)])
}
  
tewAnalyze = tew %>% 
  filter(grepl('Plas', block) | grepl('NA12878', block)) %>% # Let's only consider the NA12878 block for now
  mutate(replicate = as.integer(str_extract(block, '[0-9]$'))) %>% 
  spread(allele, count) %>% 
  group_by(snp) %>% 
  summarise
  summarize(meanAPlasmid = mean(count[grepl('Plasmid', block) & allele == 'A']),
            meanBPlasmid = mean(count[grepl('Plasmid', block) & allele == 'B']))
