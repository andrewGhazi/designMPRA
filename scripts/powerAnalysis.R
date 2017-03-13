# Let's try to redo the analysis on Ulirsch with a DESeq2 model

library(tidyverse)
library(magrittr)
library(DESeq2)

setwd('~/designMPRA/')

#### Load and reshape the counts into a matrix ####
dat = read_tsv('data/RBC_MPRA_minP_raw.txt',
               col_types = cols(chr = col_character()))

constructDat = dat %>% 
  group_by(construct) %>% 
  summarise(nRef = sum(type == 'Ref'),
            nMut = sum(type == 'Mut')) %>% 
  filter(nRef == 14, nRef == nMut)


dat %<>% filter(construct %in% constructDat$construct) # Take out the constructs that don't have 14 barcodes per

reshapeFun = function(constructDat){
  constructDat %>% 
    group_by(type) %>% 
    mutate(barcodeNum = 1:length(type)) %>% 
    ungroup %>% 
    gather(libr, count, contains('NA')) %>% 
    unite(typeLibNum, type, libr, barcodeNum) %>% 
    select(-byallele) %>% 
    spread(typeLibNum, count)
}

dat %<>% 
  group_by(construct) %>% 
  do(reshapeFun(.)) %>% 
  ungroup

mat = dat %>% 
  select(contains('K562')) %>% 
  as.data.frame %>% 
  as.matrix

colDat = mat %>% 
  colnames %>% 
  data.frame(acid = str_extract(., '[DR]NA'),
             type = str_extract(., 'Ref|Mut'),
             source = str_extract(., 'CTRL|GATA1|DNA'))
