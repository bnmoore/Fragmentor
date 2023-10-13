library(shiny)
library(data.table)
library(tidyverse)
library(Spectra)
library(rhandsontable)

shinyServer(function(input, output) {
  
  #Parse sequence string into the mono mass
  getMass = function(sequenceString, mono){
    
    sequence = grep("[A-z]+", sequenceString, value = TRUE)
    mods = grep("([0-9]+)", sequenceString, value = TRUE)
    
    chars = data.table(Abbrev1 = unlist(strsplit(sequence, split = "")))

    #A[B, on = 'a', bb := i.b]
    
    #Add up the aa letters
    if(mono){
      chars[aminoAcids, on = "Abbrev1", MonoisotopicMass := MonoisotopicMass]
      mass = sum(chars$MonoisotopicMass)
    }
    else{
      chars[aminoAcids, on = "Abbrev1", AverageMass := AverageMass]
      mass = sum(chars$AverageMass)
    }
    
    #Add in the manual mods
    
    format(mass, nsmall = 5)
  }
  
  output$elemental_comp = renderText({
    paste("Elemental Composition:", input$sequence)
  })
   
  output$theoreticalIonTable = renderDT({
    #Parse sequence
    Ntermbuild = substring(input$sequence, 1, 2:(nchar(input$sequence))-1)
    Ctermbuild = substring(input$sequence, 2:(nchar(input$sequence)), nchar(input$sequence))
    
    fragment_types = input$fragment_types
    
    #Generate ions
    #all_ions = tibble()
    #all_ions$sequence = ""
    #all_ions[atoms$Abbrev1] = 0
    
    all_ions = rbind(
      data.frame(sequence = Ntermbuild, ion_type = "a"),
      data.frame(sequence = Ntermbuild, ion_type = "b"),
      data.frame(sequence = Ntermbuild, ion_type = "c"),
      data.frame(sequence = Ctermbuild, ion_type = "x"),
      data.frame(sequence = Ctermbuild, ion_type = "y"),
      data.frame(sequence = Ctermbuild, ion_type = "z")
    )
    
    all_ions = all_ions %>% 
      left_join(ionTypeList, by = "ion_type") %>% 
      rowwise() %>% 
      mutate(ion_type = paste0(ion_type, str_length(sequence))) %>% 
      mutate(mass = getMass(sequence, input$monoOrAvg == "Monoisotopic"))
               
    #Return DT
    datatable(select(all_ions, c("ion_type", "mass")), options = list(lengthChange = FALSE))
  })
  
})






# #Parse sequence string into the mono mass
# getMass = function(sequenceString){
#   
#   sequence = grep("[A-z]+", sequenceString, value = TRUE)
#   mods = grep("([0-9]+)", sequenceString, value = TRUE)
#   
#   chars = data.table(Abbrev1 = unlist(strsplit(sequence, split = "")))
#   
#   #A[B, on = 'a', bb := i.b] 
#   chars[aminoAcids, on = "Abbrev1", MonoisotopicMass := MonoisotopicMass]
#   
#   #Add up the aa letters
#   monoMass = sum(chars$MonoisotopicMass)
#   #avgMass = sum(chars$MolecularWeight)
#   
#   #Add in the manual mods
#   
# }

# output$theoreticalIonTable = renderDT({
#   #Parse sequence
#   Ntermbuild = substring(input$sequence, 1, 2:(nchar(input$sequence))-1)
#   Ctermbuild = substring(input$sequence, 2:(nchar(input$sequence)), nchar(input$sequence))
#   
#   fragment_types = input$fragment_types
#   
#   #Generate ions
#   all_ions = rbind(
#     data.frame(sequence = Ntermbuild, ion_type = "a"),
#     data.frame(sequence = Ntermbuild, ion_type = "b"), 
#     data.frame(sequence = Ntermbuild, ion_type = "c"),
#     data.frame(sequence = Ctermbuild, ion_type = "x"),
#     data.frame(sequence = Ctermbuild, ion_type = "y"),
#     data.frame(sequence = Ctermbuild, ion_type = "z")
#   )
#   
#   all_ions = all_ions %>% 
#     left_join(ionTypeList, by = "ion_type") %>% 
#     rowwise() %>% 
#     mutate(ion_type = paste0(ion_type, str_length(sequence))) %>% 
#     mutate(mass = getMass(sequence) + 
#              MonoisotopicMass(list("C" = C,"H" = H,"N" = N,"O" = O,"S" = S,"P" = P,"Br" = Br,"Cl" = Cl,"F" = F,"Si" = Si)))
#   
#   #Return DT
#   datatable(select(all_ions, c("ion_type", "mass")), options = list(lengthChange = FALSE))
# })
# 
# })
