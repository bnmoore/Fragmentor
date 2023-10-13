library(shiny)
library(DT)
library(rhandsontable)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Fragmentor"),
  
  fluidRow(
    column(width = 1, "Nterm"),
    column(width = 8, textAreaInput("sequence", "Sequence", "RGYALG", width = "80%")),
    column(width = 1),
    column(width = 1, "Cterm"),
    column(width = 1, "Polarity")
  ),
  
  fluidRow(
    column(width = 4, textOutput("elemental_comp")),
    column(width = 4, "Deuterium Exchange")
  ),
  
  fluidRow(
    column(width = 1,
           radioButtons("monoOrAvg", label = "", choices = c("Monoisotopic", "Average")),
           checkboxGroupInput("fragment_types", "Ion types:", 
                              choiceNames = ionTypeList$ion_type, 
                              choiceValues = ionTypeList$ion_type),
           "Sidechains",
           "Losses",
           "Charge States"
           ),
    
    column(width = 8,
           DTOutput("theoreticalIonTable")
           ),
    
    column(width = 3,
           "Mass Delta area"
    )
  )
))
