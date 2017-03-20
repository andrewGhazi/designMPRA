# Let's try to redo the analysis on Ulirsch with a DESeq2 model

library(tidyverse)
library(magrittr)
library(stringr)
library(DESeq2)

setwd('~/designMPRA/')

#### Load and reshape the counts into a matrix ####
dat = read_tsv('data/RBC_MPRA_minP_raw.txt',
               col_types = cols(chr = col_character()))

rightNDat = dat %>% 
  group_by(construct) %>% 
  summarise(nRef = sum(type == 'Ref'),
            nMut = sum(type == 'Mut')) %>% 
  filter(nRef == 14, nRef == nMut)


dat %<>% filter(construct %in% rightNDat$construct) # Take out the constructs that don't have 14 barcodes per

reshapeFun = function(constructDat){
  constructDat %>% 
    group_by(type) %>% 
    mutate(barcodeNum = 1:length(type)) %>% 
    ungroup %>% 
    gather(libr, count, contains('NA')) %>% 
    unite(typeLibNum, type, libr, barcodeNum) %>% 
    select(-byallele, -clean) %>% 
    spread(typeLibNum, count)
}

dat %<>% 
  group_by(construct) %>% 
  do(reshapeFun(.)) %>% 
  ungroup

mat = dat %>% 
  select(contains('K562')) %>% 
  select(matches('CTRL|DNA')) %>% #Let's just use the control condition experiments for now
  as.data.frame %>% 
  as.matrix
rownames(mat) = dat$construct

colDat = mat %>% 
  colnames %>% 
  data_frame(column = .,
             acid = str_extract(column, '[DR]NA') %>% factor(levels = c('DNA', 'RNA')),
             type = str_extract(column, 'Ref|Mut') %>% factor(levels = c('Ref', 'Mut'))) %>% 
  as.data.frame

#### Do DESeq ####

deseqInput = DESeqDataSetFromMatrix(mat,
                                    colData = colDat %>% select(-column),
                                    design = ~acid*type)

deseqOut = DESeq(deseqInput)
rownames(deseqOut) = rownames(mat)
deseqRes = results(deseqOut)

resNames = rownames(deseqRes)
deseqRes %<>% 
  as.data.frame %>% 
  as.tbl %>% 
  mutate(construct = resNames)

#### Let's see how well these correspond to the wilcox tests ####
load('data/gatheredMPRA.RData')
goodCon = MPRA.qnactivity %>% 
  filter(!grepl('GATA', Block)) %>% 
  group_by(construct) %>% 
  summarise(good = (sum(type == 'Ref') > 0 & sum(type == 'Mut') > 0)) %>% 
  filter(good == TRUE)

MPRA.wilcox = MPRA.qnactivity %>% 
  filter(!grepl('GATA', Block)) %>% 
  filter(construct %in% goodCon$construct) %>% 
  group_by(construct) %>% 
  summarise(wilcoxP = wilcox.test(qnact[type == 'Ref'], qnact[type == 'Mut'])$p.value,
            q = p.adjust(wilcoxP, 'fdr'))

deseqRes %>% left_join(MPRA.wilcox, by = 'construct') %>% na.omit %$% cor(log(padj), log(q)) # .566537 So pretty well. Okay good.
tmp = deseqRes %>% left_join(MPRA.wilcox, by = 'construct') %>% na.omit
p = deseqRes %>% 
  left_join(MPRA.wilcox, by = 'construct') %>% 
  na.omit %>% 
  ggplot(aes(padj, q)) + 
  geom_point(alpha = .2) + 
  scale_y_log10() + 
  scale_x_log10() + 
  ggtitle('DESeq q-value vs Activity Wilcox test q-value') + 
  xlab('DESeq FDR level') + 
  ylab('Ref/Mut Activity Wilcox Test q-value')
p
ggsave(p, filename = 'outputs/QvalueComparison.png', width = 6, height = 4)

#### Let's do some bootstrapping to estimate how many barcodes are needed to detect a given log2FoldChange ####
# Let's separately quantile normalize just the ctrl block for now

MPRA.activity.mat = MPRA.activity %>% filter(construct %in% goodCon$construct, grepl('CTRL', Block)) %>% 
  select(Block, barcode, byallele, type, act) %>% 
  spread(Block, act) %>% 
  select(byallele, contains('activ')) %>% 
  as.data.frame 
rn = MPRA.activity.mat$byallele %>% as.character()
cn = MPRA.activity.mat %>% colnames %>% .[2:7]
MPRA.activity.mat %<>% .[,2:7] %>% as.matrix(dimnames = list(rn, cn))
nBarcodes = c(3,5,8,13,21,34, 40, 55, 89)






