library(shiny)
library(ratesci)


# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {


  # Compute the forumla text in a reactive expression since it is
  # shared by the output$caption and output$mpgPlot expressions
  formulaText <- reactive({
    options(digits = 4)
    out <- scoreci(x1 = as.numeric(input$x1),
                   x2 = as.numeric(input$x2),
                   n1 = as.numeric(input$n1),
                   n2 = as.numeric(input$n2),
                   dist = input$dist,
                   contrast = input$contrast,
                   skew = input$skew,
                   theta0 = theta0,
                   plot = TRUE)
    myci <- paste0("(", signif(out$estimates[1], 4), ", ",  signif(out$estimates[3], 4), ")")
    paste0("CI for ", input$x1,"/", input$n1,
           ifelse(input$contrast == "p", ":\n", paste0(" vs ", input$x2,"/", input$n2, ":\n")),
           myci)
  })

  # Return the formula text for printing as a caption
  output$caption <- renderText({
    formulaText()
  })

  # Generate a plot of the requested variable against mpg and only
  # include outliers if requested
  output$scorePlot <- renderPlot({
    scoreci(x1 = as.numeric(input$x1),
              x2 = as.numeric(input$x2),
              n1 = as.numeric(input$n1),
              n2 = as.numeric(input$n2),
              dist = input$dist,
              contrast = input$contrast,
            skew = input$skew,
            theta0 = theta0,
            plot = TRUE)
  })
})


