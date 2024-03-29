library(shiny)
library(rhandsontable)
library(shinyStorePlus)

shinyUI(fluidPage(
  
  #Initialize shinyStorePlus
  initStore(),
  
  # Application title
  fluidRow(
    column(width = 10, titlePanel("Fragmentor")),
    column(width = 2, "v2024-02-16"),
  ),
  
  fluidRow(
    column(width = 2, selectInput("Nterm_input", "Nterm:", 
                                         choices = filter(term_ref, term == "N")$terminus)),
    column(width = 6, textAreaInput("sequence_input", "Sequence", "RGYALG", 
                                    width = "100%", height = "100px")),
    column(width = 2, selectInput("Cterm_input", "Cterm:", 
                                         choices = filter(term_ref, term == "C")$terminus)),
    column(width = 2, actionButton("polarity_input", label = "+"))
  ),
  
  fluidRow(
    column(width = 8, 
           htmlOutput("elemental_comp"),
           htmlOutput("aa_comp"),
           htmlOutput("manual_mods") ),
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
           checkboxGroupInput("sidechain_input", "Sidechain losses:", 
                              choices = c("M*-sidechain", "fragment-sidechain")),
           numericInput("charge_input", label = "Charge", value = 1),
           radioButtons("max_charge_input", label = "", choices = c("Max Charge", "Selected Charge"))
           ),
    
    column(width = 6,
           rHandsontableOutput("ion_hot")
           ),
    
    column(width = 4,
           numericInput("mass_tol_input", label = "Mass Tolerance (m/z)", value = 0.4),
           radioButtons("mass_delta_input", label = "", choices = c("ppm", "m/z")),
           checkboxInput("off_by_one_input", label = "Allow +/- 1 Dalton match"),
           rHandsontableOutput("search_hot"),
           br(),
           htmlOutput("sequence_picture")
    ),
    
    #Temp data storage
    conditionalPanel('false', textInput("search_hot_store", "", ",,,,,,,,,"))
  )
))
