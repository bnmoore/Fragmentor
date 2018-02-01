library(shiny)

shinyServer(function(input, output) {
   
  output$theoreticalIonTable = renderDT(
    iris, options = list(lengthChange = FALSE)
  )
  
})
