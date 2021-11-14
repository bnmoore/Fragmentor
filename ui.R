library(shiny)
library(DT)
# library(rsconnect)
# deployApp()

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Fragmentor"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("fragmentTypes", "Ion types:", choiceNames = ionTypeList, choiceValues = ionTypeList
      )
    ),
    
    mainPanel(
      textAreaInput("sequence", "Sequence", "RGYALG", width = "80%"),
      DTOutput("theoreticalIonTable")
    )
  )
))
