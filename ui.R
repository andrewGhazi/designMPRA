
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
                                               value = 1.1),
                                  numericInput('nsnp',
                                               label = h4('Number of SNPs tested in assay'),
                                               value = 8000),
                                  numericInput('nBlock',
                                               label = h4('Number of transfections'),
                                               value = 5),
                                  numericInput('nbarcode',
                                               label = h4('Number of barcodes per SNP'),
                                               value = 14),
                                  numericInput('alpha',
                                               label = h4('Bonferroni corrected alpha'),
                                               value = .05,
                                               min = 2e-16,
                                               max = 1)),
                 conditionalPanel(condition = 'input.selectedTab==3',
                                  numericInput('nBCperSNP',
                                               label = h4('Number barcodes per SNP'),
                                               value = 14),
                                  numericInput('contextWidth',
                                               label = h4('Genomic context width (1/2 construct width)'),
                                               value = 75),
                                  textInput('fwprimer',
                                            label = h4('Forward PCR primer'),
                                            value = 'ACTGGCCAG'),#From Namrata 8/11 - original CTGGCCAGTG # old 9-21 ACTGGCCAGTG
                                  textInput('revprimer',
                                            label = h4('Reverse PCR primer'),
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
                 br(),
                 p('Here you can inspect the interaction between the various parameters of a MPRA assay and their effect on the power the detect transcriptional shifts of a given size.'),
                 h4('Inputs'),
                 tags$div(
                   tags$ul(
                     tags$li('Activity standard deviation - depending on the SNP in question and its underlying mechanism of action, typical values range from .3 to 2.'),
                     tags$li('Number of SNPs - the total number of SNPs tested in the assay.'),
                     tags$li('Number of transfections - the number of repeated transfections to act as biological replicates'),
                     tags$li('Number of barcodes per SNP - the number of barcode replicates for each SNP'),
                     tags$li('Bonferroni corrected alpha - alpha level cutoff used for calling variants as \"functional\".')
                   )
                 ),
                 p('Repeated transfections combined with quantile normalization to remove batch effects can be an effective way of improving power, however this comes at the price of increased sequencing needs. Large batch effects can reduce the statistical gains from repeated transfections, but for the purposes of this visualization they are treated as independent and identically distributed measurements.'),
                 h4('Output'),
                 p('This plot shows the power of a t-test to distinguish differential activity in a variant at a Bonferonni corrected alpha = .05 level. Note that activity (and thus transcriptional shift) are log values, and thus a transcriptional shift of +4 means the mutant allele produces exp(+5) = 148.4 times as much mRNA per input DNA molecule. While activity measurements do not strictly satisfy the assumptions of a t-test (they are not normally distributed), in practice a normal approximation usually works well.'),
                 br(),
                 plotOutput('powerPlot'),
                 value = 2),
        tabPanel("Create sequences",
                 p('Here you can upload a VCF to automatically generate the necessary construct sequences for your own MPRA experiment. The default values are comparable to those used in Ulirsch et al. 2016, but determine what works best for your needs on the Power tab.'),
                 h3('Inputs'),
                 h4('Parameters'),
                 tags$div(
                   tags$ul(
                     tags$li('Genomic context width - the amount of genomic context on either side of input SNPs.'),
                     tags$li('Number of barcodes per SNP'),
                     tags$li('PCR primers')
                 )),
                 p(),
                 h4('VCFs'),
                 p('Only the CHROM, POS, REF, and ALT columns are used. The INFO column is used only for detecting reverse strand constructs.'),
                 p('Current input constraints are: '),
                 tags$div(
                   tags$ul(
                     tags$li('Restriction sites are limited to XbaI, KpnI, and SfiI as they are used in MPRA protocol published in Melnikov (2014) (which is intended to be used with the pMPRA1 vector series).'),
                     tags$li('Barcodes are constrained to being 12bp in length'),
                     tags$li('The current total number of barcodes (the number of SNPs times the number of barcodes per SNP times 2 for reference/alternate alleles) can be at most 1,140,292.'),
                     tags$li('Insertions and deletions must encode the reference and mutant alleles (respectively) as a dash character \'-\'.'),
                     tags$li('Multiple alternate alleles should be separated in the ALT field by a comma and no spaces'),
                     tags$li('By default, the program pulls the sequence context from the forward (+) strand of the reference genome. If the user wishes to generate SNPs for genes that normally are read from the reverse strand, add a string containing "MPRAREV" to the INFO field of the VCF.')
                   )
                 ),
                 p('VCFs generated by batch querying rsID\'s on dbSNP should meet most of the formatting requirements. However the MPRAREV tag will need to be added by the user (where appropriate) because the VCF\'s do not always specify which strand the relevant gene is on.'),
                 h3('Outputs'),
                 p('Once processing has completed, click the download button to download a tab-separated file containing the sequences and other relevant information. Details:'),
                 tags$div(
                   tags$ul(
                     tags$li('sequence context comes from the hg38 assembly'),
                     tags$li('Some snps may fail due to their sequence context containing or generating an undesired restriction site')
                   )
                 ),
                 br(),
                 p('This diagram shows the layout of each MPRA construct sequence:'),
                 br(),
                 img(src = 'MPRAseqDiagram.png',
                     align = 'center',
                     width = 955,
                     height = 20),
                 br(),
                 p(''),
                 p('The outputs of this application are purely to make designing MPRA experiments more convenient for researchers. We make no guarantee as to the accuracy of the outputs and highly encourage you to methodically check your sequences before synthesizing your construct library.'),
                 br(),
                 h4('Job Time Estimate:'),
                 htmlOutput('timeText'),
                 h4('Input head:'),
                 tableOutput('inputHead'),
                 h4('Failed snps:'),
                 tableOutput('failed'),
                 value = 3),
        id = 'selectedTab'
      )
    )
    
  )

    
  )
)
