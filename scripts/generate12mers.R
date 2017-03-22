#It's gotten too tricky to randomly sample tags then throw out those that fit
#our constraints, so we're just going to start with the set of all 10bp
#sequences and cut it down, then sample from that.

# This is modified from suitableTags.R in the plateletMPRA directory

library(tcR)
library(Biostrings)
nucruns = vector(mode = 'character', length = 4) %>% DNAStringSet
ni = 1
for (i in 4) {
  for (j in c('A', 'G', 'T', 'C')) {
    nucruns[ni] = rep(j, i) %>% paste(collapse = '') %>% DNAStringSet
    ni = ni + 1
  }
}

twelvemers = generate.kmers(12) %>% DNAStringSet
cat(paste0('done generating 12mers at ', Sys.time()))

#tmp = twelvemers[sample.int(length(twelvemers), size = 10)]

#Each nucleotide occurs at least one
missingone = apply(alphabetFrequency(twelvemers)[,1:4], 1, function(x){any(x == 0)})
twelvemers = twelvemers[!missingone]
cat(paste0('done removing twelvemers missing a nucleotide at ', Sys.time()))


#Cut out those with nucleotide runs of 4 or more in a row
hasnucruns = vcountPDict(nucruns, twelvemers) %>% colSums
hasnucruns = hasnucruns > 0
twelvemers = twelvemers[!hasnucruns]
cat(paste0('done removing 12mers with runs of 4 or more at ', Sys.time()))

#Cut out those that start with TCT (creates an alternative digestion site for XbaI)
tctStart = subseq(twelvemers, 1, 3) == DNAString('TCT')
twelvemers = twelvemers[!tctStart]
cat(paste0('done removing 12mers starting with TCT at ', Sys.time()))

#Cut out those that match the miRNA seed sequences 
#For now let's just use the human ones since there are fewer and it won't take as long
source('~/plateletMPRA/mirBaseMunging.R')
haveSeedlist = vwhichPDict(humanSeedSeqs, twelvemers) #this takes ~40 minutes. All seeds takes ~1h45m
save(list = c('twelvemers', 'haveSeedlist', 'humanSeedSeqs'), file = '~/designMPRA/outputs/haveHumanRNAiSeeds.RData')
haveSeed = sapply(haveSeedlist, function(x){length(x) > 0})

twelvemers = twelvemers[!haveSeed]
cat(paste0('done removing those with mirSeeds at ', Sys.time()))

print(length(twelvemers))
save(twelvemers, file = '~/designMPRA/outputs/inertTwelveMers.RData')
