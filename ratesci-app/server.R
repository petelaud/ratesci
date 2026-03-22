library(shiny)
library(ratesci)


# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {

  formulaText1 <- reactive({
    paste0(as.numeric(input$level), "% CI for ",
           input$contrast, ": ", input$x1,"/", input$n1,
           ifelse(input$contrast == "p", ":  ",
                  paste0(" vs ", input$x2,"/", input$n2, ":  ")))
  })


  formulaText2 <- reactive({
    options(digits = 4)
    out <- scoreci(x1 = as.numeric(input$x1),
                   x2 = as.numeric(input$x2),
                   n1 = as.numeric(input$n1),
                   n2 = as.numeric(input$n2),
                   level = as.numeric(input$level) / 100,
                   dist = input$dist,
                   contrast = input$contrast,
                   skew = input$skew,
                   #                   theta0 = as.numeric(input$theta0),
                   plot = FALSE,
                   precis = as.numeric(input$precis))
    myci <- paste0("(",
                   formatC(out$estimates[1],
                           format = ifelse(input$contrast %in% c("p", "RD"), "f", "fg"),
                           as.numeric(input$precis),
                           flag = "#"),
                   ", ",
                   formatC(out$estimates[3],
                           format = ifelse(input$contrast %in% c("p", "RD"), "f", "fg"),
                           digits = as.numeric(input$precis),
                           flag = "#"),
                   ")")
    paste0(myci)
  })

  # Return the formula text for printing as a caption
  output$caption1 <- renderText({
    formulaText1()
  })

  output$caption2 <- renderText({
    formulaText2()
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
#            theta0 = as.numeric(input$theta0),
            plot = TRUE,
            precis = as.numeric(input$precis))
  })


})


