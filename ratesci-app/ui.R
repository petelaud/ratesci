library(shiny)

# Define UI for ratesci application
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Skewness-Corrected Asymptotic Score confidence intervals for comparison of rates"),

  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("dist", "Distribution:",
                list("Binomial" = "bin",
                     "Poisson" = "poi")),
    selectInput("contrast", "Contrast:",
                list("Rate Difference" = "RD",
                     "Rate Ratio" = "RR",
                     "Odds Ratio" = "OR",
                     "Proportion" = "p")),
    textInput("x1", "x1:", "10"),
    textInput("n1", "n1:", "100"),
    conditionalPanel('input.contrast != "p"',
      textInput("x2", "x2:", "10"),
      textInput("n2", "n2:", "100")
    ),
    textInput("level", "Confidence level (%):", "95"),
    selectInput("skew", "Skewness correction:", list("TRUE", "FALSE")),
#    conditionalPanel('input.contrast == "p"',
#                     selectInput("theta0", "Null hypothesis theta:", "0.5")
#    ),
#    conditionalPanel('input.contrast == "RD"',
#                     selectInput("theta0", "Null hypothesis theta:", "0")
#    ),
#    conditionalPanel('input.contrast == "RR" || input.contrast == "OR"',
#                     selectInput("theta0", "Null hypothesis theta:", "1")
#    ),
    width = 3

#    checkboxInput("outliers", "Show outliers", FALSE)
  ),

  mainPanel(
    h3(textOutput("caption")),

    plotOutput("scorePlot")
  )
))
