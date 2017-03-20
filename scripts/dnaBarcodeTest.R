# Let's try out this DNA barcodes package Ed pointed out to me

# source("https://bioconductor.org/biocLite.R")
# biocLite("DNABarcodes")

library(DNABarcodes)
library(Biostrings)
library(tidyverse)
select = dplyr::select

myTenmers = create.dnabarcodes(10)

pilotSeqs = read_csv('~/plateletMPRA/outputs/Tue_Sep_27_10-46-44_2016Pilot_81_controls.csv') %>% 
 rbind(read_csv('~/plateletMPRA/outputs/Thu_Sep_22_10-11-45_2016Pilot_81_CD36.csv') %>% select(-rs))

# analyse.barcodes(pilotSeqs$barcode) 
# > analyse.barcodes(pilotSeqs$barcode)
# Description   hamming  seqlev levenshtein
# 1               Mean Distance  7.490018 5.04016    6.305496
# 2             Median Distance  8.000000 5.00000    6.000000
# 3            Minimum Distance  1.000000 1.00000    1.000000
# 4            Maximum Distance 10.000000 8.00000   10.000000
# 5 Guaranteed Error Correction  0.000000 0.00000    0.000000
# 6  Guaranteed Error Detection  0.000000 0.00000    0.000000

# Seems okay?



# myTenmers = create.dnabarcodes(10,
#                                heuristic = 'ashlock',
#                                cores = 10)

# Great. They've put a lot more work into ensuring the robustness of a barcode
# set, so I'm just going to generate a huge set of 12-mers, filter them based on
# the constraints I was using in ~/plateletMPRA/suitableTags.R that are from
# Ulirsch (which shouldn't cut out too many) then use that.
