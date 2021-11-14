library(shiny)
library(data.table)

shinyServer(function(input, output) {
  
  #Parse sequence string into the mono mass
  getMass = function(sequenceString){
    
    sequence = grep("[A-z]+", sequenceString, value = TRUE)
    mods = grep("([0-9]+)", sequenceString, value = TRUE)
    
    chars = data.table(Abbrev1 = unlist(strsplit(sequence, split = "")))
    
    #A[B, on = 'a', bb := i.b] 
    chars[aminoAcids, on = "Abbrev1", MonoisotopicMass := MonoisotopicMass]
    
    #Add up the aa letters
    monoMass = sum(chars$MonoisotopicMass)
    #avgMass = sum(chars$MolecularWeight)
    
    #Add in the manual mods
    
  }
   
  output$theoreticalIonTable = renderDT({
    #Parse sequence
    Ntermbuild = substring(input$sequence, 1, 2:(nchar(input$sequence))-1)
    Ctermbuild = substring(input$sequence, 2:(nchar(input$sequence)), nchar(input$sequence))
    
    input$fragmentTypes
    
    #Generate ions
    
    #Return datatable
    all_ions = rbind(data.frame(sequence = Ntermbuild), data.frame(sequence = Ctermbuild))
    all_ions = all_ions %>% 
      rowwise() %>% 
      mutate(mass = getMass(sequence))
    
    datatable(all_ions, options = list(lengthChange = FALSE))
  })
  
})
