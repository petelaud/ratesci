library(shiny)
library(ratesci)


# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {

  formulaText <- reactive({
    options(digits = 4)
    out <- scoreci(x1 = as.numeric(input$x1),
                   x2 = as.numeric(input$x2),
                   n1 = as.numeric(input$n1),
                   n2 = as.numeric(input$n2),
                   level = as.numeric(input$level) / 100,
                   dist = input$dist,
                   contrast = input$contrast,
                   skew = input$skew,
                   theta0 = as.numeric(input$theta0),
                   plot = TRUE)
    myci <- paste0("(", signif(out$estimates[1], 4), ", ",
                   signif(out$estimates[3], 4), ")")
    paste0(as.numeric(input$level), "% CI for ", input$contrast, ": ", input$x1,"/", input$n1,
           ifelse(input$contrast == "p", ":  ",
                  paste0(" vs ", input$x2,"/", input$n2, ":  ")),
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
            level = as.numeric(input$level) / 100,
            dist = input$dist,
              contrast = input$contrast,
            skew = input$skew,
            theta0 = as.numeric(input$theta0),
            plot = TRUE)
  })
})


