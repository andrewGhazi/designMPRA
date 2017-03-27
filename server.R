
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(tidyverse)
library(stringr)
library(pwr)
#library(VariantAnnotation)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(magrittr)
#source('~/designMPRA/scripts/processVCF.R')
source('~/designMPRA/scripts/processVCFfast.R')
#expand = S4Vectors::expand

load('~/designMPRA/outputs/inertTwelveMers.RData')
mers = twelvemers

shinyServer(function(input, output) {

  output$powerPlot <- renderPlot({

    plotDat = data_frame(meanDiff = seq(0,5, by = .05),
                         pwr = pwr.t.test(n = input$nbarcode*input$nBlock,
                                          d = meanDiff / input$sigma,
                                          sig.level = input$alpha / input$nsnp)$power) # This is where the bonferonni correction somes in
    
    ggplot(plotDat, aes(meanDiff, pwr)) + 
      geom_line() + 
      xlab('Variant Transcriptional Shift (absolute value)') + 
      ylab('Power to detect difference with given design') + 
      ylim(0,1) + 
      theme(text = element_text(size = 14))

  })
  
  inVCF = reactive({
    inVCF = input$vcf
    
    if (is.null(inVCF)) {
      return(NULL)
    }
    
    # Screw it, we're going back to using readr. VariantAnnotation and GRanges are too irritating to work with.
    skipNum = system(paste0('grep ^## ', inVCF$datapath, ' | wc -l'), intern = TRUE) %>% as.numeric
    
    #Check that the header doesn't have spaces in place of tabs. If it does, replace the spaces with tabs and create a new col_names variable
    vcfColumns = system(paste0('head -', skipNum + 1, ' ', inVCF$datapath, ' | tail -1'), 
                        intern = TRUE) %>% 
      gsub('#', '', .) %>% 
      gsub('[ ]+', '\t', .) %>% #replace spaces with tabs if applicable
      str_split('\t') %>% 
      unlist
    
    read_tsv(inVCF$datapath, 
             skip = skipNum + 1,
             col_names = vcfColumns)
  })
  
  output$testHead = eventReactive(input$Go, {
    head(mtcars)
  })
  
  output$timeText = renderText({
    if (is.null(inVCF())) {
      return(NULL)
    }
    
    nSnp = inVCF() %>% nrow
    
    timeGuess = round(nSnp*input$nBCperSNP*10/1000, digits = 3)
    nSnpStatement = paste0('You are requesting ', 
                           input$nBCperSNP, 
                           ' barcoded sequences for each of ', 
                           nSnp, 
                           ' snps. This yields a total of ', 
                           nSnp*input$nBCperSNP, ' 
                             sequences. Each sequence takes roughly 10ms to be generated, so this request will require roughly <b>', 
                           timeGuess,
                           ' seconds </b>to process.')
    
    nSnpStatement
  })
  
  vcfOut = eventReactive(input$Go, {
    print('We\'re in!')
    
    validate(
      need(input$nBCperSNP > 0, 
           'Please input a number of barcodes to use per SNP')
      )
    
    validate(
      need(input$nBCperSNP*nrow(inVCF())*2 < 20000,
           'Your request will take too long (see the estimated time above). Try running the package version of this software locally.')
    )
    
    style = 'notification'
    progress = shiny::Progress$new(style = style, min = 0, max = 1)
    progress$set(message = 'Generating sequences...', value = 0)
    
    
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    res = processVCF(inVCF(), input$nBCperSNP, input$contextWidth, input$fwprimer, input$revprimer, updateProgress)
    progress$close
    return(res)
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
      head
  })
  
  output$downloadSequences = downloadHandler(
    filename = function() {'out.tsv'},
    content = function(file){
      write_tsv(vcfOut()$result, file)
    }
  )

})
