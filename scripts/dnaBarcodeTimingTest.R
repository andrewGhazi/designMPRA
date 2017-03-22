# Creating DNA barcodes was pretty slow for 10mers, so let's look at how the ashlock timing increases as you go from just 4 to 7


library(DNABarcodes)
library(Biostrings)
library(tidyverse)
select = dplyr::select

barcodeFun = function(n){
  #Creates the barcodes and returns a list containing them & how long it took
  strt = Sys.time()
  
  bcs = create.dnabarcodes(n,
                           heuristic = 'conway', #ran this the first time with ashlock
                           cores = 10)
  
  stp = Sys.time()
  bcTime = stp - strt
  
  res = list(bcTime = bcTime, barcodes = bcs)
  return(res)
}

times = data_frame(n = 4:12,
                   barcodeSet = n %>% map(barcodeFun),
                   time = barcodeSet %>% map_dbl(~.x$bcTime),
                   setSize = barcodeSet %>% map_int(~length(.x$barcodes)))
save(times, file = 'outputs/dnaBarcodeTimingTestConway.RData')

p = ggplot(times, aes(n, time)) + 
  geom_line() +
  scale_y_log10()
ggsave('outputs/dnaBarcodeTimingConway.png')

# > times
# # A tibble: 9 × 5
# n barcodeSet         time setSize      seconds
# <int>     <list>        <dbl>   <int>        <dbl>
#   1     4 <list [2]> 7.779598e-04      12 7.779598e-04
# 2     5 <list [2]> 9.925365e-04      30 9.925365e-04
# 3     6 <list [2]> 2.973795e-03      76 2.973795e-03
# 4     7 <list [2]> 3.806424e-02     265 3.806424e-02
# 5     8 <list [2]> 2.017884e-01     701 2.017884e-01
# 6     9 <list [2]> 4.330076e+00    2595 4.330076e+00
# 7    10 <list [2]> 5.633678e+01    9154 5.633678e+01
# 8    11 <list [2]> 7.547426e+00   26702 4.528440e+02
# 9    12 <list [2]> 2.513768e+00   98527 9.049320e+03
# The set sizes are too small :/
