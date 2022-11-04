#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(gghighlight)
library(ggplot2)
library(ggvis)
library(ggsci)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  load("human_data.RData")
  
  obsv_diso <- reactive({paste(input$predictor,"_psites_obsv_diso",sep = "")})
  expct_diso <- reactive({paste(input$predictor,"_psites_expct_diso",sep = "")})
  binom_sig <- reactive({paste(input$predictor,"_binom_sig",sep = "")})
  padj <- reactive({paste(input$predictor,"_binom_q",sep = "")})
  
  kinase_column <-reactive({input$kinase_column})
  
  kinase_target <- reactive({
    
    if (input$kinase_column=="target") {"Cdk1 target"}
    else if (input$kinase_column=="target_mapk") {"mapk target"}
    else if (input$kinase_column=="target_aurk") {"aurk target"}
    else if (input$kinase_column=="target_plk") {"plk target"}
    else if (input$kinase_column=="target_nek") {"nek target"}
    else if (input$kinase_column=="target_dyrk") {"dyrk target"}
    
  })
  
  
  output$kinasesTable <- renderDataTable({
    
    
    if (input$kinase_column!="NULL") {
      aux_table <- subset(human_data,get(kinase_column()) == kinase_target())
      aux_table <- aux_table[,c("ACC#","GENE","PROTEIN", "length",expct_diso(),obsv_diso(),padj())]
      colnames(aux_table) <- c("UniProt ID","Gene","Protein","Length","Expected Psites in IDRs","Observed Psites in IDRs","Adjusted P-value")
      aux_table
    }
    else {
      aux_table <- human_data[,c("ACC#","GENE","PROTEIN", "length",expct_diso(),obsv_diso(),padj())]
      colnames(aux_table) <- c("UniProt ID","Gene","Protein","Length","Expected Psites in IDRs","Observed Psites in IDRs","Adjusted P-value")
      aux_table
    }},
    
    options = list(
      pageLength = 10,
      FixedHeader = TRUE,
      search = F,
      tabIndex= F
      
    ))
  
  
  output$distPlot <- renderPlot({
    all_phospho_plot<-ggplot(human_data) +
      geom_point(aes(x=.data[[obsv_diso()]],y=.data[[expct_diso()]], colour = .data[[binom_sig()]]),size=2,alpha=0.80)+
      geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
      ggpubr::theme_classic2() +
      theme(text = element_text(size=15),legend.position = c(0.25,0.75),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
      guides(color=guide_legend(title="Statistical significance")) +
      scale_x_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ xlab("Observed phospho S/T in IDR") +
      scale_y_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ ylab("Expected phospho S/T in IDR") +
      scale_colour_manual(values = c(pal_jco()(10)[3],"#ffdd15ff",pal_jco()(10)[4]))
    
    # .data[[]] functions converts from the string to the variable name in the inherited data frame. In this case human_data.
    if (input$kinase_column!="NULL") {all_phospho_plot + gghighlight(.data[[kinase_column()]] == kinase_target(),keep_scales =T ,unhighlighted_params = list(color="grey",alpha=0.3))}
    else {all_phospho_plot}
  },width = 400,height = 400)
  
}




#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



# Define UI for application that draws a histogram
ui <- fluidPage(
  
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
  
  
  
)


shinyApp(ui = ui, server = server)
