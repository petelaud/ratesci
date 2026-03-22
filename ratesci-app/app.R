
library(shiny)
library(ratesci)

# Define UI for ratesci application
ui <- pageWithSidebar(

  # Application title
  headerPanel("Skewness-Corrected Asymptotic Score confidence intervals for comparison of rates"),

  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("dist", "Distribution:",
                list("Binomial" = "bin",
                     "Poisson" = "poi")),
    conditionalPanel('input.dist == "bin"',
                     selectInput("contrast", "Contrast:",
                                 list(
                                   "Rate Difference" = "RD",
                                   "Rate Ratio" = "RR",
                                   "Odds Ratio" = "OR",
                                   "Proportion" = "p")),
    ),
    conditionalPanel('input.dist == "poi"',
                     selectInput("contrast", "Contrast:",
                                 list("Rate Difference" = "RD",
                                      "Rate Ratio" = "RR",
                                      "Proportion" = "p"))
    ),
    textInput("x1", "x1:", "10"),
    textInput("n1", "n1:", "100"),
    conditionalPanel('input.contrast != "p"',
                     textInput("x2", "x2:", "10"),
                     textInput("n2", "n2:", "100"),
    ),
    textInput("level", "Confidence level (%):", "95"),
    selectInput("skew", "Skewness correction:", list("TRUE", "FALSE")),
    textInput("precis", "Decimal precision:", "4"),
    width = 3

  ),

  mainPanel(
    h3(textOutput("caption1")),
    h3(textOutput("caption2")),
    plotOutput("scorePlot")
  )
)


#    conditionalPanel('input.contrast == "p"',
#                     selectInput("theta0", "Null hypothesis theta:", "0.5")
#    ),
#    conditionalPanel('input.contrast == "RD"',
#                     selectInput("theta0", "Null hypothesis theta:", "0")
#    ),
#    conditionalPanel('input.contrast == "RR" || input.contrast == "OR"',
#                     selectInput("theta0", "Null hypothesis theta:", "1")
#    ),



server <- function(input, output) {

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


}

shinyApp(ui = ui, server = server)
