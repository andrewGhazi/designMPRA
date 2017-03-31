#Intended to be run from ~/designMPRA/
filesToDeploy = c('server.R',
                  'ui.R',
                  'scripts/processVCFfast.R',
                  'outputs/inertTwelveMersChar.RData',
                  'www/formulas.png',
                  'www/MPRAdiagram.png',
                  'www/MPRAseqDiagram.png')

deployApp(appFiles = filesToDeploy)