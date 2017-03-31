
library(shiny)
library(tidyverse)
library(stringr)
library(pwr)
#library(VariantAnnotation)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(magrittr)

dir = '~/designMPRA/data/CD36_final86_dbSNP.vcf'

skipNum = system(paste0('grep ^## ', dir, ' | wc -l'), intern = TRUE) %>% as.numeric

#Check that the header doesn't have spaces in place of tabs. If it does, replace the spaces with tabs and create a new col_names variable
vcfColumns = system(paste0('head -', skipNum + 1, ' ', dir, ' | tail -1'), 
                    intern = TRUE) %>% 
  gsub('#', '', .) %>% 
  gsub('[ ]+', '\t', .) %>% #replace spaces with tabs if applicable
  str_split('\t') %>% 
  unlist

vcf = read_tsv(dir, 
         skip = skipNum + 1,
         col_names = vcfColumns)

nper = 5
seqwidth = 75
fwprimer = 'ACTGGCCAG'
revprimer = 'CTCGGCGGCC'
