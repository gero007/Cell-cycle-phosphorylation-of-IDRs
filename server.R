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
library(ggvis)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    load("plotting_data_test.RData")
    output$distPlot <- renderPlot({

      obsv_diso <- paste(input$predictor,"_psites_obsv_diso",sep = "")
      expct_diso <- paste(input$predictor,"_psites_expct_diso",sep = "")
      binom_sig <- paste(input$predictor,"_binom_sig",sep = "")
        
      all_phospho_plot<-ggplot(human_data) +
        geom_point(aes(x=.data[[obsv_diso]],y=.data[[expct_diso]], colour = .data[[binom_sig]]),size=2,alpha=0.80)+
        geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
        ggpubr::theme_classic2() +
        theme(text = element_text(size=15),legend.position = c(0.25,0.75),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
        guides(color=guide_legend(title="Statistical significance")) +
        scale_x_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ xlab("Observed phospho S/T in IDR") +
        scale_y_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ ylab("Expected phospho S/T in IDR") +
        scale_colour_manual(values = c(pal_jco()(10)[3],"#ffdd15ff",pal_jco()(10)[4]))
        
        # .data[[]] functions converts from the string to the variable name in the inherited data frame. In this case human_data.
        if (input$kinase_column!="NULL") {all_phospho_plot + gghighlight(.data[[input$kinase_column]] == "Cdk1 target",keep_scales =T ,unhighlighted_params = list(color="grey",alpha=0.3),)} else {all_phospho_plot}
        },width = 400,height = 400)
        
})
