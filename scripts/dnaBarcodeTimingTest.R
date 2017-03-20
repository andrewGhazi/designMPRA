# Creating DNA barcodes was pretty slow for 10mers, so let's look at how the ashlock timing increases as you go from just 4 to 7


library(DNABarcodes)
library(Biostrings)
library(tidyverse)
select = dplyr::select

barcodeFun = function(n){
  #Creates the barcodes and returns a list containing them & how long it took
  strt = Sys.time()
  
  bcs = create.dnabarcodes(n,
                           heuristic = 'ashlock',
                           cores = 10)
  
  stp = Sys.time()
  bcTime = stp - strt
  
  res = list(bcTime = bcTime, barcodes = bcs)
  return(res)
}

times = data_frame(n = 4:9,
                   barcodeSet = n %>% map(barcodeFun),
                   time = barcodeSet %>% map_dbl(~.x$bcTime))
save(times, file = 'outputs/dnaBarcodeTimingTest.RData')

p = ggplot(times, aes(n, time)) + 
  geom_line()
ggsave('outputs/dnaBarcodeTiming.png')