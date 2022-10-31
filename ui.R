#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Cell cycle phosphorylation of IDRs"),
    sidebarPanel( 
    
      selectInput("predictor", "Disorder predictor: ",
                  c("SPOT disorder" = "SPOT",
                    "IUPred" = "IUPred"
                  )),
      
      selectInput("kinase_column", "Kinase:",
                c("All phosphorylations" = "NULL",
                  "CDK1 subfamily" = "target",
                  "MAPK family" = "target_mapk",
                  "AURK family" = "target_aurk",
                  "PLK family" = "target_plk",
                  "NEK family" = "target_nek",
                  "DYRK family" = "target_dyrk"
                  )),
      
      width = 2),
    

    # Show a plot of the generated distribution
    mainPanel(dataTableOutput("kinasesTable"),
              plotOutput("distPlot")
              )
    
    
    
))
