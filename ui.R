library(shiny)
library(DT)
# library(rsconnect)
# deployApp()

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Fragmentor"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("fragmentTypes", "Choose fragment types:",
                         choiceNames =
                           list("a", "b", "c", "x", "y", "z"),
                         choiceValues =
                           list("a", "b", "c", "x", "y", "z")
      )
    ),
    
    mainPanel(
      textAreaInput("sequence", "Sequence", "RGYALG", width = "1000px"),
      DTOutput("theoreticalIonTable")
    )
  )
))
