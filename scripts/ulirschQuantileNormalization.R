#Let's redo the first part's of Ulirsch's analysis so we can get a reportable number on the activity std dev

library(tidyverse)
library(magrittr)
library(stringr)
library(DESeq2)
select = dplyr::select

setwd('~/designMPRA/')

#### Load and reshape the counts into a matrix ####
dat = read_tsv('data/RBC_MPRA_minP_raw.txt',
               col_types = cols(chr = col_character()))


dat %<>% 
  dmap_at(grep('K562', names(dat)), ~ .x*1e6/sum(.x)) %>% # counts per million
  mutate(DNA = (K562_minP_DNA1 + K562_minP_DNA2)/2) %>% #take mean of DNA's
  select(-contains('GATA1'), -contains('minP_DNA')) %>% 
  filter(DNA != 0)

RNA = dat %>% select(contains('RNA'))
activities = RNA %>% dmap(~log(.x / dat$DNA))

qnact = activities %>% as.data.frame() %>% as.matrix %>% normalize.quantiles()

names(activities) %<>% gsub('minP_', '', .)
groupSD = dat %>% select(-contains('K562'), -DNA) %>% 
  bind_cols(activities) %>% 
  gather(block, value, contains('K562')) %>% 
  filter(is.finite(value)) %>% 
  group_by(block, byallele) %>% 
  summarise(stdDev = sd(value))

groupStdDevLimits = groupSD %>% 
  group_by(block) %>% 
  summarise(SDecdf = ecdf(stdDev) %>% list,
            lowerBound = SDecdf %>% map_dbl(~quantile(.x, .05)),
            med = SDecdf %>% map_dbl(~quantile(.x, .5)),
            upperBound = SDecdf %>% map_dbl(~quantile(.x, .95)),
            foldChange = exp(med)) # This is an important result to keep at hand
save(groupStdDevLimits, file = 'outputs/groupStdDevLimits.RData')

ggplot(groupSD, aes(stdDev)) + 
  geom_histogram(bins = 100) + 
  facet_grid(block ~ .)  + # So these are the distributions of activities in the different blocks
  ggtitle('Activity standard deviations by block in Ulirsch et al., 2016') + 
  theme(strip.text = element_text(size = 6))
ggsave('outputs/UlirschStdDevByBlock.png')

#You can see that the activities are typically around 1 within a single block

activities %<>% gather(block, activity)

activities %>% 
  ggplot(aes(block, activity)) + 
  geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

activities %>% filter(!is.infinite(activity)) %>%  group_by(block) %>% summarise(stdDev = sd(activity))

# So the spread in the quantile normalized control groups is much tighter

#Looking at the result from Ulirsch's original analysis:
load("~/MPRA/origWorkspace.RData")
MPRA_minP.melt %>% as.tbl() %>% filter(grepl('CTRL', variable)) %>% group_by(byallele) %>% summarise(stdDev = sd(value)) %>% ggplot(aes(stdDev)) + geom_density()
MPRA_minP.melt %>% as.tbl() %>% filter(grepl('CTRL', variable)) %>% group_by(byallele) %>% summarise(stdDev = sd(value)) %>% .$stdDev %>% mean
