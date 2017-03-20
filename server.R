
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(tidyverse)
library(stringr)
library(pwr)
library(VariantAnnotation)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
source('scripts/processVCF.R')


shinyServer(function(input, output) {

  output$powerPlot <- renderPlot({

    plotDat = data_frame(meanDiff = seq(0,5, by = .05),
                         pwr = pwr.t.test(n = input$nbarcode,
                                          d = meanDiff / input$sigma,
                                          sig.level = .05 / input$nsnp)$power) # This is where the bonferonni correction somes in
    
    ggplot(plotDat, aes(meanDiff, pwr)) + 
      geom_line() + 
      xlab('Transcriptional Shift (absolute value)') + 
      ylab('Power to detect difference in 1 transfection') + 
      ylim(0,1)

  })
  
  inVCF = reactive({
    inVCF = input$vcf
    
    if (is.null(inVCF)){
      return(NULL)
    }
    
    # skipNum = system(paste0('grep ^## ', inVCF$datapath, ' | wc -l'), intern = TRUE) %>% as.numeric
    # 
    # #Check that the header doesn't have spaces in place of tabs. If it does, replace the spaces with tabs and create a new col_names variable
    # vcfColumns = system(paste0('head -', skipNum + 1, ' ~/plateletMPRA/data/CD36_initial75_dbSNP.vcf | tail -1'), 
    #                     intern = TRUE) %>% 
    #   gsub('#', '', .) %>% 
    #   gsub('[ ]+', '\t', .) %>% #replace spaces with tabs if applicable
    #   str_split('\t') %>% 
    #   unlist
    
    readVcf(inVCF$datapath, 'hg38')
    # read_tsv(inVCF$datapath, 
    #          skip = skipNum + 1,
    #          col_names = vcfColumns)
  })
  
  vcfOut = eventReactive(input$Go, {
    inVCF() %>% processVCF(input$nBCperSNP, input$contextWidth)
  })

  
  output$inputHead = renderTable({
    if (is.null(inVCF())){
      return(NULL)
    }
    
    inVCF() %>% 
      rowRanges %>% 
      as.data.frame %>% 
      mutate(.,
             ALT = ALT %>% map_chr(~as.character(paste(.x, collapse = ','))), 
             snp = rownames(.)) %>% 
      head
  })
  
  output$outputHead = eventReactive(input$Go,{
    validate(
      need(input$nBCperSNP > 0,
           'Please input a number of barcodes per snp to proceed')
    )
    
    validate(
      need(nrow(inVCF()) > 0,
           'Please input a VCF')
    )
    
    validate(
      need(nrow(vcfOut) > 0,
           'Sequence generation failed')
    )
    
    head(vcfOut())
  })
  
  output$downloadSequences = downloadHandler(
    'out.tsv',
    content = function(file){
      write.csv(vcfOut(), file)
    })

})
