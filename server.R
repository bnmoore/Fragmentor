library(shiny)
library(data.table)
library(tidyverse)
library(rhandsontable)
#library(rawrr)

shinyServer(function(input, output) {
  
  # RGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALG
  
  
#Functions
  get_mass = function(table){
    if(input$mono_input == "Monoisotopic")
      table$MonoisotopicMass
    else
      table$AverageMass
  }
  

  
  
#Reactive functions  
  all_ions = reactive({
    #if(is.null(input$fragment_types_input))
      #return(data.frame())
    
    
    sequence_input = input$sequence_input
    
    #Parse sequence
    Mbuild = sequence_input
    Ntermbuild = substring(sequence_input, 1, 2:(nchar(sequence_input))-1)
    Ctermbuild = substring(sequence_input, 2:(nchar(sequence_input)), nchar(sequence_input))
    
    selected_ion_types = ion_types_ref[0,]
    
    for(i in input$fragment_types_input){
      selected_ion_types = rbind(selected_ion_types, filter(ion_types_ref, ion_type == i))
    }
    
    selected_ion_types_M = filter(selected_ion_types, term == "M")
    selected_ion_types_N = filter(selected_ion_types, term == "N")
    selected_ion_types_C = filter(selected_ion_types, term == "C")
    
    #Generate ions
    all_ions = rbind(
      expand.grid(sequence = Mbuild, 
                  ion_type = selected_ion_types_M$ion_type, 
                  charge = seq(1, 3, 1), 
                  terminus = "M", 
                  stringsAsFactors = FALSE),
      expand.grid(sequence = Ntermbuild, 
                  ion_type = selected_ion_types_N$ion_type, 
                  charge = seq(1, 3, 1), 
                  terminus = input$Nterm_input, 
                  stringsAsFactors = FALSE),
      expand.grid(sequence = Ctermbuild, 
                  ion_type = selected_ion_types_C$ion_type, 
                  charge = seq(1, 3, 1), 
                  terminus = input$Cterm_input, 
                  stringsAsFactors = FALSE)
    )
    
    #Sequence + Iontype + Terminus + Mods + Charge --> Atoms --> Mass
    a0 = data.frame(t(sapply(unique(all_ions$sequence), sequence_to_atoms)), check.names = FALSE) %>% 
      mutate_all(as.numeric)
    a0$sequence = rownames(a0)
    all_ions0 = all_ions %>% 
      left_join(a0, by = "sequence")
    
    all_ions1 = all_ions %>% 
      left_join(ion_types_ref, by = "ion_type")
    
    all_ions2 = all_ions %>% 
      filter(terminus != "M") %>% 
      left_join(term_ref, by = "terminus") 

    all_ions3 = all_ions %>% 
      filter(terminus == "M") %>% 
      mutate(terminus = input$Nterm_input)%>% 
      left_join(term_ref, by = "terminus") 
    
    all_ions4 = all_ions %>% 
      filter(terminus == "M") %>% 
      mutate(terminus = input$Cterm_input)%>% 
      left_join(term_ref, by = "terminus") 

    all_ions5 = all_ions %>%
      mutate(`H+` = charge)
    
    all_ions98 = bind_rows(all_ions0, all_ions1, all_ions2, all_ions3, all_ions4, all_ions5)
    
    all_ions = all_ions98 %>%
      ungroup() %>% 
      group_by(sequence, ion_type, charge) %>% 
      summarise(across(C:`e-`, ~ sum(.x, na.rm = TRUE)))
    
    polarity = input$polarity_input
    
    all_ions = all_ions %>% 
      ungroup() %>% 
      rowwise() %>% 
      mutate(ion_name = paste0(ion_type, "<sub>", str_length(sequence), "</sub> ", charge, polarity)) %>% 
      mutate(ion_type_charge = paste0(ion_type, " ", charge, polarity)) %>% 
      mutate(length = str_length(sequence))
    
    #Mass calculation
    mass = all_ions %>% 
      select(C:`e-`)
    mass = as.matrix(mass) %*% diag(atom_table_wide)
    mass = rowSums(mass)
    all_ions = cbind(all_ions, mass)
    all_ions = all_ions %>% 
      mutate(mass = mass / charge) %>% 
      arrange(mass) %>% 
      mutate(mass = format(mass, nsmall = 5))
    
    
    if(input$table_view_input == "List"){
      all_ions = all_ions %>% 
        select(c("ion_name", "mass"))
    }
    else{
      all_ions = all_ions %>% 
        select(c("ion_type_charge", "mass", "length")) %>% 
        arrange(ion_type_charge) %>% 
        pivot_wider(names_from = length, values_from = mass)
    } 
    all_ions
  })

  
#Main code
  output$elemental_comp = renderText({
    paste("Elemental Composition:", " ")
  })
  
  output$atom_comp = renderText({
    paste("Atom Composition:", " ")
  })
   
  #Change this to reactive on sequence input
  output$ion_table = renderRHandsontable({
    rhandsontable(all_ions(), rowHeaders = FALSE, readOnly = TRUE) %>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
      hot_cols(colWidths = 100) %>% 
      hot_cols(renderer = htmlwidgets::JS("safeHtmlRenderer")) 
  })
  
  output$mass_table = renderRHandsontable({
    mass_table_df = data.frame(`Mass` = 0, `Intensity` = 0, `Int Error` = 0, `Ion` = "", `Mass Delta` = 0, check.names = FALSE)
    
    rhandsontable(mass_table_df, rowHeaders = FALSE)%>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
  })
  
})









# sequence_to_mass = function(sequenceString){
#   sequenceString = grep("[A-z]+", sequenceString, value = TRUE)
#   mods = grep("([0-9]+)", sequenceString, value = TRUE)
#   
#   chars = data.table(Abbrev1 = unlist(strsplit(sequenceString, split = "")))
#   
#   #Add up the aa letters
#   if(input$mono_input == "MonoisotopicMass"){
#     chars[data.table(amino_acids_ref), on = "Abbrev1", MonoisotopicMass := MonoisotopicMass]
#     mass = sum(chars$MonoisotopicMass)
#   }
#   else{
#     chars[data.table(amino_acids_ref), on = "Abbrev1", AverageMass := AverageMass]
#     mass = sum(chars$AverageMass)
#   }
#   
#   #Add in the manual mods
# }
# charge_carrier_mass = function(){
#   if(input$polarity_input == "+"){
#     out = atoms_ref %>% filter(Abbrev1 == "H+") %>% 
#       select(MonoisotopicMass, AverageMass)
#   }
#   if(input$polarity_input == "-"){
#     out = atoms_ref %>% filter(Abbrev1 == "H+") %>% 
#       select(MonoisotopicMass, AverageMass) %>% 
#       mutate(MonoisotopicMass = -MonoisotopicMass) %>% 
#       mutate(AverageMass = -AverageMass)
#   }
#   get_mass(out)
# }



