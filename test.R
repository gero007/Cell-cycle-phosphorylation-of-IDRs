library(shiny)
library(plotly)
library(dplyr)

ui <- fluidPage(
  actionButton("go", "Click to recalc"),
  plotlyOutput("plot")
)

gendata <- function(){
  ndata <- 10
  d <- tibble(text=LETTERS[1:ndata], f=1, x=runif(ndata)) %>% mutate(r = rank(x))
  rbind(mutate(d, x=-1), d, mutate(d, x=-1)) %>%
    arrange(text)
}

server <- function(input, output, session){
  
  origdata <- gendata()
  if (FALSE){ # for offline testing
    print(head(origdata))
    my <- list(olddata = origdata, newdata = origdata)
  }
  
  my <- reactiveValues(
    olddata = origdata,
    newdata = origdata
  )
  
  output$plot <- renderPlotly({
    cat("renderPlotly\n")
    plot_ly() %>%
      add_trace(x=origdata$x, y=origdata$r, frame=origdata$f, line=list(width=20, simplify=FALSE), type="scatter", opacity=0.5, mode="lines", name="Rank") %>%
      add_trace(x=origdata$x + 0.02, y=origdata$r, frame=origdata$f, text=origdata$text, type="scatter", mode="text", showlegend=FALSE) %>%
      layout(xaxis=list(range=list(0,1.1))) %>%
      animation_opts(frame=500, transition=500, redraw=FALSE)
  })
  
  observeEvent(input$go, {
    req(my$newdata)
    cat("observeEvent input$go\n")
    my$olddata <- my$newdata # save old data
    my$newdata <- gendata() %>% # generate new data
      mutate(f=my$olddata$f+1)
    print(head(my$newdata))
    # https://plot.ly/javascript/plotlyjs-function-reference/#plotlyanimate
    plotlyProxy("plot", session=session, deferUntilFlush=FALSE) %>%
      plotlyProxyInvoke("animate",
                        # frameOrGroupNameOrFrameList
                        list(
                          data = list(list(
                            x = my$newdata$x,
                            y = my$newdata$r,
                            frame = my$newdata$f
                          ),
                          list(
                            x = my$newdata$x + 0.02,
                            y = my$newdata$r,
                            text = my$newdata$text,
                            frame = my$newdata$f
                          )),
                          traces = list(as.integer(0), as.integer(1)),
                          layout = list()
                        ),
                        # animationAttributes
                        list()
      )# plotlyProxyInvoke
  })
  
}

shinyApp(ui, server)