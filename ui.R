
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("MPRA design tools"),

  sidebarLayout(
    sidebarPanel(conditionalPanel(condition = 'input.selectedTab==2',
                                  numericInput('sigma',
                                               label = h4('Activity standard deviation'),
                                               value = 1.6),
                                  numericInput('nsnp',
                                               label = h4('Number of snps tested in assay'),
                                               value = 8000),
                                  numericInput('nbarcode',
                                               label = h4('n barcodes per snp'),
                                               value = 20)),
                 conditionalPanel(condition = 'input.selectedTab==3',
                                  numericInput('nBCperSNP',
                                               label = h6('n barcodes per snp'),
                                               value = 0),
                                  numericInput('contextWidth',
                                               label = h6('Genomic context width (1/2 construct width)'),
                                               value = 75),
                                  textInput('fwprimer',
                                            label = h6('Forward PCR primer'),
                                            value = 'ACTGGCCAG'),#From Namrata 8/11 - original CTGGCCAGTG # old 9-21 ACTGGCCAGTG
                                  textInput('revprimer',
                                            label = h6('Reverse PCR primer'),
                                            value = 'CTCGGCGGCC'),
                                  fileInput('vcf',
                                            label = 'Input VCF'),
                                  actionButton('Go',
                                               label = 'Go!'),
                                  downloadButton('downloadSequences', 'Download Output'))),
    mainPanel(
      tabsetPanel(
        tabPanel("MPRA", 
                 p("Massively parallel reporter assays (MPRA) are high-throughput assays designed to inspect the transcriptional effects of many variants in parallel. See the diagram below for a graphical explanation."),
                 br(),
                 p("The \"activity\" of an allele is defined as the log-ratio of the output mRNA count to the input DNA count for a given allele. Because there is significant noise in these counts, allele constructs are replicated across multiple barcodes. The observed activity of each barcode is taken as one observation of the activity of the allele. The \"transcriptional shift\" of a variant is defined as the difference in activity between the mutant and reference alleles. Thus variants that increase transcription (e.g. the red variant in the diagram below) have a positive transcriptional shift and those that decrease transcription (e.g. the blue variant) have a negative shift."),
                 br(),
                 img(src = 'formulas.png',
                     align = 'center',
                     width = 728,
                     height = 102),
                 p("This application aims to help researchers understand the relationships between power, noise, and experimental design in MPRAs. Optimizing experimental design is critical for sensitivity in a MPRA assay -- the variance of how many transcripts a single cell makes from a single DNA molecule can be extremely high. Multiple testing corrections further reduce sensitivity. Repeated transfections combined with quantile normalization to remove batch effects can be an effective way of improving power, however this comes at the price of increased sequencing needs."),
                 img(src = 'MPRAdiagram.png',
                     align = 'center',
                     width = 723,
                     height = 947),
                 p(''),
                 br(),
                 p(''),
                 value = 1),
        tabPanel("Power", 
                 p('Here you can inspect the interaction between the various parameters of a MPRA assay and their effect on the power the detect transcriptional shifts of a given size. Note that activity (and thus transcriptional shift) are log values, and thus a transcriptional shift of +4 means the mutant allele produces exp(+5) = 148.4 times as much mRNA per input DNA molecule.'),
                 br(),
                 p('Repeated transfections combined with quantile normalization to remove batch effects can be an effective way of improving power, however this comes at the price of increased sequencing needs.'),
                 br(),
                 p('This plot shows the power of a t-test to distinguish differential activity in a variant at a Bonferonni corrected alpha = .05 level. While activity measurements do not strictly satisfy the assumptions of a t-test (they are not normally distributed), in practice a normal approximation usually works well.'),
                 br(),
                 plotOutput('powerPlot'),
                 value = 2),
        tabPanel("Create sequences",
                 p('Here you can upload a VCF to automatically generate the necessary construct sequences for your own MPRA experiment. Only the CHROM, POS, REF, and ALT columns are used. Note that the input sequence context width is HALF of the total context (e.g. input of 75 yields 150bp of sequence context). '),
                 br(), 
                 p('Current constraints are: '),
                 tags$div(
                   tags$ul(
                     tags$li('sequence context comes from the hg38 assembly'),
                     tags$li('Restriction sites are limited to XbaI, KpnI, and SfiI as they are used in MPRA protocol published in Melnikov (2014) (which is intended to be used with the pMPRA1 vector series).'),
                     tags$li('Some snps may fail due to their sequence context containing or generating an undesired restriction site'),
                     tags$li('Barcodes are constrained to being 10bp in length'),
                     tags$li('The current total number of barcodes can be at most 131587. This is currently being addressed to increase this number to ~16 million.'),
                     tags$li('Insertions and deletions must encode the reference and mutant alleles (respectively) as a dash character \'-\'.')
                   )
                 ),
                 p('The outputs of this application are purely to make designing MPRA experiments more convenient for researchers. We make no guarantee as to the accuracy of the outputs and highly encourage you to methodically check your sequences before synthesizing your construct library.'),
                 br(),
                 h4('Input head:'),
                 tableOutput('inputHead'),
                 tableOutput('testHead'),
                 value = 3),
        id = 'selectedTab'
      )
    )
    
  )

    
  )
)
