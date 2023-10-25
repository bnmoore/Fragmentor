library(shiny)
library(rhandsontable)

shinyUI(fluidPage(
  
  # Application title
  fluidRow(
    column(width = 10, titlePanel("Fragmentor")),
    column(width = 2, "v2023-10-25")
  ),
  
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
    column(width = 8, 
           htmlOutput("elemental_comp"),
           htmlOutput("aa_comp") ),
    column(width = 4, checkboxInput("deuterium_exchange_input", label = "Deuterium Exchange"))
  ),
  
  br(),
  
  fluidRow(
    column(width = 2,
           radioButtons("mono_input", label = "", choices = c("Monoisotopic", "Average")),
           radioButtons("table_view_input", label = "", choices = c("List", "Table")),
           checkboxGroupInput("fragment_types_input", "Ion types:", 
                              choices = filter(ion_types_ref, !str_detect(ion_type, "'"))$ion_type, 
                              selected = c("M","b","y")),
           checkboxGroupInput("losses_input", "Losses:", 
                              choices = unique(losses_ref$loss)),
           numericInput("charge_input", label = "Charge", value = 1),
           radioButtons("max_charge_input", label = "", choices = c("Max Charge", "Selected Charge"))
           ),
    
    column(width = 6,
           rHandsontableOutput("ion_hot")
           ),
    
    column(width = 4,
           rHandsontableOutput("search_hot")
    )
  )
))
