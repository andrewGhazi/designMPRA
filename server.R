
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

load('~/designMPRA/outputs/inertTwelveMers.RData')
mers = twelvemers

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
  
  output$timeText = renderText({
    if (is.null(inVCF())) {
      return(NULL)
    }
    
    nSnp = inVCF() %>% expand %>% nrow
    
    timeGuess = round(nSnp*input$nBCperSNP*80/1000, digits = 3)
    nSnpStatement = paste0('You are requesting ', 
                           input$nBCperSNP, 
                           ' barcoded sequences for each of ', 
                           nSnp, 
                           ' snps. This yields a total of ', 
                           nSnp*input$nBCperSNP, ' 
                             sequences. Each sequence takes roughly 80ms to be generated, so this request will require roughly <b>', 
                           timeGuess,
                           ' seconds </b>to process.')
    
    
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

  output$failed = renderTable({
    data.frame(snpID = vcfOut()$failed)
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
