library(shiny)
library(rhandsontable)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Fragmentor"),
  
  fluidRow(
    column(width = 2, selectInput("Nterm_input", "Nterm:", 
                                         choices = filter(term_ref, term == "N")$terminus)),
    column(width = 6, textAreaInput("sequence_input", "Sequence", "RGYALG", 
                                    width = "100%", height = "100px")),
    column(width = 2, selectInput("Cterm_input", "Cterm:", 
                                         choices = filter(term_ref, term == "C")$terminus)),
    column(width = 2, radioButtons("polarity_input", label = "", choices = c("+", "-")))
  ),
  
  fluidRow(
    column(width = 4, textOutput("elemental_comp")),
    column(width = 4, textOutput("atom_comp")),
    column(width = 4, checkboxInput("deuterium_exchange_input", label = "Deuterium Exchange"))
  ),
  
  br(),
  
  fluidRow(
    column(width = 2,
           radioButtons("mono_input", label = "", choices = c("Monoisotopic", "Average")),
           radioButtons("table_view_input", label = "", choices = c("List", "Table")),
           checkboxGroupInput("fragment_types_input", "Ion types:", 
                              choices = ion_types_ref$ion_type, 
                              selected = c("M", "b","y")),
           "Sidechains",
           "Losses",
           "Charge States"
           ),
    
    column(width = 6,
           rHandsontableOutput("ion_table")
           ),
    
    column(width = 4,
           rHandsontableOutput("mass_table")
    )
  )
))
