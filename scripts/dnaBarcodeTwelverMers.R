# It will take roughly 21 hours to generate the 12mers on 10 cores. Let's do it, and with 12 cores.

library(DNABarcodes)
library(Biostrings)
library(tidyverse)
select = dplyr::select

myTwelveMers = create.dnabarcodes(12,
                               heuristic = 'ashlock',
                               cores = 12,
                               iterations = 110)

save(myTwelveMers, file = 'outputs/12mersRaw.RData')