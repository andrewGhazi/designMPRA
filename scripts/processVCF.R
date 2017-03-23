#Function to take an input vcf and generate an output vcf with the specified context width and 

processVCF = function(vcf, nper, seqwidth, fwprimer, revprimer){
  library(BSgenome.Hsapiens.UCSC.hg38)
  
  expand = S4Vectors::expand
  select = dplyr::select
  
  load('~/designMPRA/outputs/inertTwelveMers.RData')
  mers = twelvemers
  
  genome = BSgenome.Hsapiens.UCSC.hg38
  info(vcf)$rs = rownames(info(vcf))
  snps = vcf %>% expand
  
  genome = BSgenome.Hsapiens.UCSC.hg38
  
  #seqwidthOG = seqwidth
  fwprimer %<>% DNAString 
  revprimer %<>% DNAString
  
  kpn = DNAString('GGTACC') #KpnI
  xba = DNAString('TCTAGA') #XbaI
  sfi = DNAString('GGCCNNNNNGGCC') #SfiI
  
  maxbc = length(snps)*2*nper
  bcsavailable = 1:length(mers)
  
  bc = vector(mode = 'character', length = maxbc) %>% DNAStringSet
  
  bc.ind = 1
  
  snplist = vector('list', maxbc) #This will contain all of the final information for each oligo
  i.list = vector('list', 7) #The info for each construct gets put into this list, then this gets assigned to the right element of snplist
  names(i.list) = c('snpdat', 'type', 'constructindex', 'totindex', 'barcode', 'seq', 'rs')
  
  for (i in 1:length(snps)) {
    
    # get +/- seqwidth position
    rangestart = start(ranges(rowRanges(snps)))[i] - seqwidth
    rangeend = end(ranges(rowRanges(snps)))[i] + seqwidth
    
    snpseq = subseq(genome[[paste0('chr', as.character(seqnames(rowRanges(snps))[i]))]], # the chrom field needs to be only digits
                    start = rangestart, 
                    end = rangeend)
    snpdat = snps %>% .[i] %>% rowRanges
    
    
    for (j in 1:2) { #1 = Ref, 2 = Mut
      
      #This if else block generates the genomic sequence with or without the mutant allele
      if (j == 1) { #if you're doing the reference sequence
        #mid = ref(snps)[[i]]
        constrseq = snpseq
        typeval = 'Ref'
      } else {#if you're doing the mutant sequence
        
        mid = alt(snps)[[i]] #get the mutant allele #I spread out the row with two alternate alleles across two rows earlier
        
        refwidth = snps %>% ref %>% width %>% .[i] #get the width of the reference allele
        refall = ref(snps)[[i]] #get the reference allele
        
        isSNV = (refall  %>% toString %in% c('A', 'G', 'T', 'C')) & (mid %>% toString %in% c('A', 'G', 'T', 'C'))
        isdel = mid %>% toString == '-' #THIS DOESN'T WORK FOR SNPS IN GENERAL, JUST THE ONES LEN SENT US
        isins = refall %>% toString == '-' #THIS DOESN'T WORK FOR SNPS IN GENERAL, JUST THE ONES LEN SENT US
        
        if (isSNV) { #if the SNP is an SNV
          constrseq = c(subseq(snpseq, 1, seqwidth), 
                        mid, 
                        subseq(snpseq, seqwidth + 2, nchar(snpseq))) #construct sequence
        } else if (isdel) { #if the SNP is a deletion
          constrseq = c(subseq(snpseq, 1, seqwidth), 
                        subseq(snpseq, seqwidth + refwidth + 1, nchar(snpseq))) #construct sequence
        } else if (isins) {  #if the SNP is an insertion
          constrseq = c(subseq(snpseq, 1, seqwidth), 
                        mid, 
                        subseq(snpseq, seqwidth + 1, nchar(snpseq))) #construct sequence
        } else{
          stop(paste0(info(expand(vcf))$rs[i], ' was not able to be classified as a SNV, deletion, or insertion'))
        }
        
        typeval = 'Mut'
      }
      
      # This block generates the specified number of constructs per allele
      for (k in 1:nper) {
        
        #Randomly choose from the barcodes available
        newbc.i = sample.int(bcsavailable %>% length, 1) 
        newbc = mers[[bcsavailable[newbc.i]]]
        
        final = c(fwprimer, #This is the statements that concatenates all of the elements together
                  DNAString('TG'), # Not sure what the purpose of this is but Len said to add it
                  constrseq, 
                  kpn, 
                  xba, 
                  newbc, 
                  DNAString('GGC'),
                  revprimer)
        
        #count how many digestion sites occur. Should only be the two we intentionally add in.
        ndigsite = sum(countPattern(kpn, final), # This should be 1
                       countPattern(xba, final), # This should be 1
                       countPattern(sfi, final, fixed = FALSE), #fixed = FALSE allows the ambiguity code in sfi to match any letter instead of just 'N' # This should be 0, also in promoter #because of the added GGC this is now 1
                       countPattern(kpn %>% rev, final), #This should be 0
                       countPattern(xba %>% rev, final), #This should be 0
                       countPattern(sfi %>% rev, final, fixed = FALSE)) # This should be 0
        
        bcattempts = 1
        while(ndigsite>3 && bcattempts < 6){ # if any extra digestion sites are created, keep trying until you only have the two you need
          #This could happen by weirdness at the junctions that is hard to account for so just keep trying until it works
          #There have to be the two we put in deliberately
          newbc.i = sample.int(bcsavailable %>% length, 1)
          newbc = mers[[bcsavailable[newbc.i]]]
          final = c(fwprimer, #This is the statements that concatenates all of the elements together
                    DNAString('TG'), # Not sure what the purpose of this is but Len said to add it
                    constrseq, 
                    kpn, 
                    xba, 
                    newbc, 
                    DNAString('GGC'),
                    revprimer)
          #lens = c(fwprimer %>% length, constrseq %>% length, kpn %>% length, xba %>% length, newbc %>% length, revprimer %>% length)
          ndigsite = sum(countPattern(kpn, final), # This should be 1
                         countPattern(xba, final), # This should be 1
                         countPattern(sfi, final, fixed = FALSE), #fixed = FALSE allows the ambiguity code in sfi to match any letter instead of just 'N' # This should be 0
                         countPattern(kpn %>% rev, final), #This should be 0
                         countPattern(xba %>% rev, final), #This should be 0
                         countPattern(sfi %>% rev, final, fixed = FALSE)) # This should be 0
          bcattempts = bcattempts + 1
          #if(bcattempts>5){stop('5 barcodes were attempted without success. Something fishy\'s goin\' on around heah\'')}
        }
        
        if (bcattempts > 5){
          
          # the snp is bad. Increment the barcode index (leaving the intermediate elements of the final snp list blank) and break
          bc.ind = bc.ind + nper
          break
          
        } else {
          # the snp is good
          
          #and at this point we know the tag used is good, so we remove it from the list of barcodes available
          bcsavailable = bcsavailable[-which(1:length(bcsavailable) == newbc.i)]
          bc[[bc.ind]] = newbc
          
          i.list[1] = snpdat #Add all the information into a list #
          i.list[2] = typeval
          i.list[3] = k
          i.list[4] = bc.ind
          i.list[5] = newbc
          i.list[6] = final
          i.list[7] = info(snps)$rs[i]
          snplist[bc.ind] = list(i.list) #and put it into our final results list
          bc.ind = bc.ind + 1
        }
        
      }
    }
    cat(i)
    #setTxtProgressBar(pb, i)
  }; #close(pb)
  
  
  #Some munging to account for snps we couldn't get because they either missed one allele or both entirely
  resDat = data_frame(startList = snplist,
                      nll = startList %>% map_lgl(~!is.null(.x))) %>% 
    filter(nll) %>% 
    mutate(rs = startList %>% map_chr(~.x$rs),
           type = startList %>% map_chr(~.x$type)) %>% 
    group_by(rs) %>% 
    summarise(bothTypes = length(type))
  
  gotBoth = resDat %>% filter(bothTypes > nper) %>% .$rs
  
  snplist = snplist[-which(sapply(snplist, is.null))]
  snplist = snplist[(snplist %>% map_chr(~.x$rs) %>% unique) %in% gotBoth]
  
  for (x in 1:length(snplist)) {
    snplist[[x]]$totindex = x
  }
  
  resdf = data_frame(res = snplist,
                     totIndex = res %>% map_int(~.x$totindex),
                     constructIndex = res %>% map_int(~.x$constructindex),
                     rs = res %>% map_chr(~.x$rs),
                     type = res %>% map_chr(~.x$type),
                     ref = res %>% map_chr(~toString(.x$snpdat$REF)),
                     alt = res %>% map_chr(~toString(.x$snpdat$ALT)),
                     barcode = res %>% map_chr(~toString(.x$barcode)),
                     seq = res %>% map_chr(~toString(.x$seq))) %>% 
    select(-res)
  
  inputRS = rownames(info(vcf))
  failedSnps = inputRS[!(inputRS %in% resdf$rs)]
  res = list(result = resdf, failed = failedSnps)
  write_tsv(resdf, paste0('~/designMPRA/outputs/seqFileOutputs/', Sys.Date() %>% gsub('-', '_', .), '.tsv'))
  return(res)
}