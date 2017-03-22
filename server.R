
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
#library(BSgenome.Hsapiens.UCSC.hg38)
library(magrittr)
source('~/designMPRA/scripts/processVCF.R')
expand = S4Vectors::expand

load('data/inertTenmers.RData') # just using these until 12mers are ready
mers = tenmers

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
    
    if (is.null(inVCF)) {
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
  })
  
  output$testHead = eventReactive(input$Go, {
    head(mtcars)
  })
  
  vcfOut = eventReactive(input$Go, {
    print('We\'re in!')
    
    validate(
      need(input$nBCperSNP > 0, 
           'Please input a number of barcodes to use per SNP')
      )
    processVCF(inVCF(), input$nBCperSNP, input$contextWidth, input$fwprimer, input$revprimer)
  })
  
  tmp = observeEvent(input$Go, {
    print('test')
    print(str(vcfOut()))
  })

  
  output$inputHead = renderTable({
    if (is.null(inVCF())) {
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
  
  output$downloadSequences = downloadHandler(
    filename = function() {'out.tsv'},
    content = function(file){
      write_tsv(vcfOut()$result, file)
    }
  )

})
