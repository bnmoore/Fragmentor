library(shiny)
library(data.table)
library(tidyverse)
library(rhandsontable)

shinyServer(function(input, output) {
  
  # RGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALG
  # ACDEFGHIKLMNPQRSsTtUVWYy

  
#Reactive functions 
  all_ions = reactive({
    if(input$sequence_input == "" | length(input$fragment_types_input) == 0)
      return(data.frame())
    
    #Gather inputs
    sequence_input = input$sequence_input
    polarity = input$polarity_input
    mono = input$mono_input == "Monoisotopic"
    deut = input$deuterium_exchange_input
    
    
    if(input$max_charge_input == "Max Charge")
      chargebuild = seq(1, input$charge_input, 1)
    else
      chargebuild = input$charge_input
      
    if(polarity == "-")
      chargebuild = -chargebuild
    
    
    #Parse sequence
    Mbuild = sequence_input
    Ntermbuild = substring(sequence_input, 1, 2:(nchar(sequence_input))-1)
    Ctermbuild = substring(sequence_input, 2:(nchar(sequence_input)), nchar(sequence_input))
    
    selected_ion_types = ion_types_ref[0,]
    
    for(i in input$fragment_types_input){
      selected_ion_types = rbind(selected_ion_types, filter(ion_types_ref, ion_type == i))
      if(i == "d")
        selected_ion_types = rbind(selected_ion_types, filter(ion_types_ref, ion_type == i) %>% mutate(ion_type = "d'"))
      if(i == "w")
        selected_ion_types = rbind(selected_ion_types, filter(ion_types_ref, ion_type == i) %>% mutate(ion_type = "w'"))
    }
    
    selected_ion_types_M = filter(selected_ion_types, term == "M")
    selected_ion_types_N = filter(selected_ion_types, term == "N")
    selected_ion_types_C = filter(selected_ion_types, term == "C")
    
    #Sidechain loss
    selected_losses = c("", input$losses_input)
    if(sum(selected_losses == "sidechain") > 0){
      selected_losses = selected_losses[selected_losses != "sidechain"]
      selected_losses = c(selected_losses, filter(losses_ref, loss == "sidechain")$loss_display)
    }
    
    #Generate ions
    all_ions = rbind(
      expand.grid(sequence = Mbuild, 
                  ion_type = selected_ion_types_M$ion_type, 
                  charge = chargebuild, 
                  terminus = "M", 
                  loss = selected_losses,
                  stringsAsFactors = FALSE),
      expand.grid(sequence = Ntermbuild, 
                  ion_type = selected_ion_types_N$ion_type, 
                  charge = chargebuild, 
                  terminus = input$Nterm_input, 
                  loss = selected_losses,
                  stringsAsFactors = FALSE),
      expand.grid(sequence = Ctermbuild, 
                  ion_type = selected_ion_types_C$ion_type, 
                  charge = chargebuild, 
                  terminus = input$Cterm_input, 
                  loss = selected_losses,
                  stringsAsFactors = FALSE)
    )
    
    all_ions = all_ions %>%
      rowwise() %>% 
      mutate(sidechain = substring(loss, nchar(loss), nchar(loss))) %>% 
      mutate(sidechain = if_else(sidechain == "", " ", sidechain)) %>% 
      filter(sidechain == " " | (ion_type != "M" & str_detect(sequence, sidechain)))
    
    #Sequence + Iontype + Terminus + Mods + Sidechains + Charge --> Atoms --> Mass
    
    #Sequences
    a0 = data.frame(t(sapply(unique(all_ions$sequence), sequence_to_atoms, deut_exchange = deut)), check.names = FALSE) %>% 
      mutate_all(as.numeric)
    a0$sequence = rownames(a0)
    all_ions0 = all_ions %>% 
      left_join(a0, by = "sequence")
    
    
    #abcxyz
    all_ions1 = all_ions %>% 
      left_join(ion_types_ref, by = "ion_type")
    
    
    #terminus
    all_ions2 = all_ions %>% 
      filter(terminus != "M") %>% 
      left_join(term_ref, by = "terminus")

    
    #M
    all_ions3 = all_ions %>% 
      filter(terminus == "M") %>% 
      mutate(terminus = input$Nterm_input)%>% 
      left_join(term_ref, by = "terminus") 
    all_ions4 = all_ions %>% 
      filter(terminus == "M") %>% 
      mutate(terminus = input$Cterm_input)%>% 
      left_join(term_ref, by = "terminus") 

    
    #dvw    
    partial_sidechains = filter(losses_ref, d == 1) %>% select(-loss) 
    partial_sidechains_prime = filter(losses_ref, d_prime == 1) %>% select(-loss) 
    full_sidechains = filter(losses_ref, v == 1) %>% select(-loss) 
    all_ions5 = all_ions %>% 
      filter(terminus != "M") %>% 
      left_join(term_ref, by = "terminus") %>% 
      mutate(last_aa = if_else(term == "N", substring(sequence, nchar(sequence), nchar(sequence)), "")) %>% 
      mutate(last_aa = if_else(term == "C", substring(sequence, 1, 1), last_aa)) %>% 
      select(-(C:D)) %>% 
      mutate(loss = "")
    all_ions51 = all_ions5 %>% 
      filter(ion_type == "d") %>% 
      inner_join(partial_sidechains, 
                by = c("last_aa" = "Abbrev1"), relationship = "many-to-many")
    all_ions51prime = all_ions5 %>% 
      filter(ion_type == "d'") %>% 
      inner_join(partial_sidechains_prime, 
                 by = c("last_aa" = "Abbrev1"), relationship = "many-to-many")
    all_ions52 = all_ions5 %>% 
      filter(ion_type == "v") %>% 
      inner_join(full_sidechains, 
                 by = c("last_aa" = "Abbrev1"), relationship = "many-to-many")
    all_ions53 = all_ions5 %>% 
      filter(ion_type == "w") %>% 
      inner_join(partial_sidechains, 
                 by = c("last_aa" = "Abbrev1"), relationship = "many-to-many")
    all_ions53prime = all_ions5 %>% 
      filter(ion_type == "w'") %>% 
      inner_join(partial_sidechains_prime, 
                 by = c("last_aa" = "Abbrev1"), relationship = "many-to-many")
    all_ions5 = bind_rows(all_ions51, all_ions51prime, all_ions52, all_ions53, all_ions53prime)
    
    
    #Other losses
    all_ions6 = all_ions %>% 
      left_join(losses_ref, by = c("loss" = "loss_display"), relationship = "many-to-many") %>% 
      group_by(sequence, ion_type, charge, loss, loss_mass) %>% 
      filter(row_number() == 1)
    
    
    #charge
    all_ions10 = all_ions %>%
      mutate(`H+` = charge) 
    
    
    all_ions98 = bind_rows(all_ions0, all_ions1, all_ions2, all_ions3, all_ions4, all_ions5, all_ions6, all_ions10)
    
    all_ions = all_ions98 %>%
      ungroup() %>% 
      group_by(sequence, ion_type, charge, loss) %>% 
      summarise(across(C:D, ~ sum(.x, na.rm = TRUE))) %>% 
      rowwise() %>% 
      filter(!(ion_type == "d" & sum(str_detect(all_ions51$sequence, paste0("^",sequence,"$"))) == 0)) %>% 
      filter(!(ion_type == "d'" & sum(str_detect(all_ions51prime$sequence, paste0("^",sequence,"$"))) == 0)) %>% 
      filter(!(ion_type == "v" & sum(str_detect(all_ions52$sequence, paste0("^",sequence,"$"))) == 0)) %>% 
      filter(!(ion_type == "w" & sum(str_detect(all_ions53$sequence, paste0("^",sequence,"$"))) == 0)) %>% 
      filter(!(ion_type == "w'" & sum(str_detect(all_ions53prime$sequence, paste0("^",sequence,"$"))) == 0))
    

    #String formatting
    all_ions = all_ions %>% 
      ungroup() %>% 
      rowwise() %>% 
      mutate(charge_display = if_else(abs(charge)>1, as.character(abs(charge)), "")) %>% 
      mutate(ion_name = 
            paste0(ion_type, "<sub>", str_length(sequence), "</sub>", loss, "<sup>", charge_display, polarity, "</sup>")) %>% 
      mutate(ion_name = if_else(ion_type == "M", 
            paste0("[M", polarity, charge_display, "H]", loss, "<sup>", charge_display, polarity, "</sup>"), 
            ion_name)) %>%
      mutate(ion_name = if_else(ion_type == "M*", 
            paste0("[M*", polarity, charge_display, "H]", loss, "<sup>", charge_display, polarity, "</sup>"), 
            ion_name)) %>%
      mutate(ion_type_charge = paste0(ion_type, loss, "<sup>", charge_display, polarity, "</sup>")) %>% 
      mutate(length = str_length(sequence))
    
    
    #Mass calculation
    atom_table_wide = if_else(mono, atom_table_wide_mono, atom_table_wide_avg)
    
    mass = all_ions %>% 
      select(C:D)
    mass = as.matrix(mass) %*% diag(atom_table_wide)
    mass = rowSums(mass)
    all_ions = cbind(all_ions, mass)
    all_ions = all_ions %>% 
      mutate(mass = mass / abs(charge)) %>% 
      arrange(mass) %>% 
      mutate(mass = format(mass, nsmall = 5))
    
    
    
    if(input$table_view_input == "List"){
      all_ions = all_ions %>% 
        select(c("mass", "ion_name"))
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
    sequence = input$sequence_input
    deut = input$deuterium_exchange_input
    
    atoms = t(data.frame(sequence_to_atoms(sequence, deut)))
    
    label = ""
    for(a in 1:ncol(atoms)){
      if(atoms[,a] > 0)
        label = paste0(label, colnames(atoms)[a], "<sub>", atoms[,a], "</sub> ")
    }
    
    paste("Elemental Composition:", label)
  })
  
  output$aa_comp = renderText({
    sequence = input$sequence_input
    tab = strsplit(sequence, split = "") %>% 
      unlist() %>% 
      table() %>% 
      data.frame()
    
    label = ""
    for(t in 1:nrow(tab)){
      label = paste0(label, tab$.[t], tab$Freq[t], " ")
    }
    
    paste("Amino Acid Composition:", label)
  })
   
  #Change this to reactive on sequence input
  output$ion_table = renderRHandsontable({
    rhandsontable(all_ions(), rowHeaders = FALSE, readOnly = TRUE) %>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
      hot_cols(colWidths = 150) %>% 
      hot_cols(renderer = htmlwidgets::JS("Handsontable.renderers.HtmlRenderer"))
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



