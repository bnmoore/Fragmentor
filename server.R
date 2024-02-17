library(shiny)
library(data.table)
library(tidyverse)
library(rhandsontable)
library(stringr)
library(shinyStorePlus)

shinyServer(function(input, output, session) {
  
  # RGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALG
  # ACDEFGHIKL-12.0MNPQRSs-3TtUVWYy+15.
  # R-1GYALG-128.1536GGG+15.

  decimal_pattern = "[+|-]\\d*\\.?\\d*"
  bad_characters = paste0(c("[^", amino_acids_ref$Abbrev1, "]"), collapse = "")
  
#Reactive functions 
  sequence_input = reactive({
    seq = str_remove_all(input$sequence_input, decimal_pattern)
    if(str_detect(seq, bad_characters))
      ""
    else
      seq
  })
  
  mods = reactive({
    sequence = input$sequence_input
    
    mods = data.frame(mod_position = integer(), mod_mass = double())
    
    matches = unlist(str_extract_all(sequence, decimal_pattern))
    for(i in 1:length(matches))
    {
      if(length(matches) > 0){
        mods = rbind(mods, data.frame(mod_position = str_locate(sequence, decimal_pattern)[1,"start"] - 1, 
                                      mod_mass = as.numeric(matches[i])))
        sequence = str_remove(sequence, decimal_pattern)
      }
    }

    mods
  })
  
  all_ions = reactive({
    req(col_highlight)
    req(row_highlight)
    
    full_sequence = sequence_input()
    
    if(full_sequence == "" | length(input$fragment_types_input) == 0)
      return(data.frame(ion_name = character(), mass = double(), term = character(), position = integer()))
    
    #Gather inputs
    mods = mods()
    polarity = POLARITY()
    mono = input$mono_input == "Monoisotopic"
    deut = input$deuterium_exchange_input
    
    if(input$max_charge_input == "Max Charge")
      chargebuild = seq(1, input$charge_input, 1)
    else
      chargebuild = input$charge_input
      
    if(polarity == "-")
      chargebuild = -chargebuild
    
    
    #Parse sequence
    Mbuild = full_sequence
    Ntermbuild = substring(full_sequence, 1, 2:(nchar(full_sequence))-1)
    Ctermbuild = substring(full_sequence, (nchar(full_sequence)):2, nchar(full_sequence))
    
    
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
    
    #Neutral losses
    selected_losses = paste0(c("", input$losses_input), " ")
    
    
    #Sidechain loss
    if(!is.null(input$sidechain_input)){
      selected_losses = c(selected_losses, filter(losses_ref, sidechain == 1)$loss_display)
    }

    # selected_losses = paste0(c("", input$losses_input), " ")
    # if(sum(selected_losses == "sidechain ") > 0){
    #   selected_losses = selected_losses[selected_losses != "sidechain "] #This space is important
    #   selected_losses = c(selected_losses, filter(losses_ref, loss == "sidechain")$loss_display)
    # }
    
    #Generate ions
    all_ions = rbind(
      expand.grid(sequence = Mbuild, 
                  ion_type = selected_ion_types_M$ion_type, 
                  charge = chargebuild, 
                  terminus = "M",
                  term = "M",
                  loss = selected_losses,
                  stringsAsFactors = FALSE),
      expand.grid(sequence = Ntermbuild, 
                  ion_type = selected_ion_types_N$ion_type, 
                  charge = chargebuild, 
                  terminus = input$Nterm_input,
                  term = "N",
                  loss = selected_losses,
                  stringsAsFactors = FALSE),
      expand.grid(sequence = Ctermbuild, 
                  ion_type = selected_ion_types_C$ion_type, 
                  charge = chargebuild, 
                  terminus = input$Cterm_input, 
                  term = "C",
                  loss = selected_losses,
                  stringsAsFactors = FALSE)
    )
    
    all_ions = all_ions %>%
      rowwise() %>% 
      mutate(sidechain = substring(loss, nchar(loss), nchar(loss))) %>% 
      mutate(sidechain = if_else(sidechain == "", " ", sidechain)) %>% 
      filter(sidechain == " " | (ion_type != "M" & str_detect(sequence, sidechain))) %>% 
      mutate(loss = str_remove_all(loss, " "))
    
    if(!("M*-sidechain" %in% input$sidechain_input)){
      all_ions = all_ions %>% 
        filter(sidechain == " " | (ion_type != "M*" & str_detect(sequence, sidechain)))
    }
    
    if(!("fragment-sidechain" %in% input$sidechain_input)){
      all_ions = all_ions %>% 
        filter(sidechain == " " | (ion_type == "M*" & str_detect(sequence, sidechain)))
    }
    
    #Sequence + Iontype + Terminus + Mods + Sidechains + Charge --> Atoms --> Mass
    
    #Sequences
    a0 = data.frame(t(sapply(unique(all_ions$sequence), sequence_to_atoms)), check.names = FALSE) %>% 
      mutate_all(as.numeric)
    a0$sequence = rownames(a0)
    all_ions0 = all_ions %>% 
      left_join(a0, by = "sequence")
    
    
    #abcxyz
    all_ions1 = all_ions %>% 
      left_join(ion_types_ref, by = c("ion_type", "term"))
    
    
    #terminus
    all_ions2 = all_ions %>% 
      filter(terminus != "M") %>% 
      left_join(term_ref, by = c("terminus", "term"))

    
    #M
    all_ions3 = all_ions %>% 
      filter(terminus == "M") %>% 
      mutate(terminus = input$Nterm_input)%>% 
      mutate(term = "N") %>% 
      left_join(term_ref, by = c("terminus", "term")) %>% 
      mutate(term = "M")
    all_ions4 = all_ions %>% 
      filter(terminus == "M") %>% 
      mutate(terminus = input$Cterm_input)%>% 
      mutate(term = "C") %>% 
      left_join(term_ref, by = c("terminus", "term")) %>% 
      mutate(term = "M")

    
    #dvw    
    partial_sidechains = filter(losses_ref, d == 1) %>% select(-loss) 
    partial_sidechains_prime = filter(losses_ref, d_prime == 1) %>% select(-loss) 
    full_sidechains = filter(losses_ref, v == 1) %>% select(-loss) 
    all_ions5 = all_ions %>% 
      filter(sidechain == " ") %>% 
      filter(term != "M") %>% 
      left_join(term_ref, by = c("terminus", "term")) %>% 
      mutate(last_aa = if_else(term == "N", substring(sequence, nchar(sequence), nchar(sequence)), "")) %>% 
      mutate(last_aa = if_else(term == "C", substring(sequence, 1, 1), last_aa)) %>% 
      select(-(C:D)) #%>% 
      #mutate(loss = "")
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
    if(deut){
      all_ions10 = all_ions %>%
        mutate(`D+` = charge) 
    }
    else{
      all_ions10 = all_ions %>%
        mutate(`H+` = charge)
    }
    
    
    all_ions98 = bind_rows(all_ions0, all_ions1, all_ions2, all_ions3, all_ions4, all_ions5, all_ions6, all_ions10)
    
    all_ions = all_ions98 %>%
      ungroup() %>% 
      group_by(sequence, ion_type, charge, loss, term) %>%     #selection is here
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
      mutate(length = str_length(sequence)) %>% 
      mutate(position = 0) %>% 
      mutate(position = ifelse(term == "N", length, position)) %>% 
      mutate(position = ifelse(term == "C", str_length(full_sequence) - length, position))
    
    
    #Mass calculation
    atom_table_wide = if_else(mono, atom_table_wide_mono, atom_table_wide_avg)
    
    mass = all_ions %>% 
      select(C:D)
    mass = as.matrix(mass) %*% diag(atom_table_wide)
    mass = rowSums(mass)
    all_ions = cbind(all_ions, mass)
    
    #Manual mods
    if(nrow(mods) > 0)
    {
      for(i in 1:nrow(mods)){
        all_ions = all_ions %>% 
          mutate(mass = ifelse((term == "N" & mods$mod_position[i] <= position) | 
                               (term == "C" & mods$mod_position[i] > position) | 
                                term == "M",
                                mass + mods$mod_mass[i], mass))
      }
    }
    
    all_ions = all_ions %>% 
      mutate(mass = mass / abs(charge)) %>% 
      arrange(mass)
    
    all_ions = all_ions %>% 
      mutate(list_row = row_number() - 1) %>% 
      mutate(list_col = 0) %>% 
      mutate(ion_type_charge = as.factor(ion_type_charge)) %>% 
      mutate(grid_row = as.numeric(ion_type_charge) - 1) %>% 
      mutate(grid_col = length)

    all_ions
  })

  col_highlight = reactiveVal(value = c())
  row_highlight = reactiveVal(value = c())
  
  cell_renderer = "function(instance, td, row, col, prop, value, cellProperties) {
    Handsontable.renderers.HtmlRenderer.apply(this, arguments);
  
    if(instance.params){
      hcols = instance.params.col_highlight
      hcols = hcols instanceof Array ? hcols : [hcols] 
      hrows = instance.params.row_highlight
      hrows = hrows instanceof Array ? hrows : [hrows] 
      
      for (i = 0; i < hcols.length; i++) { 
        if (hcols[i] == col && hrows[i] == row) {
            td.style.background = 'lightgreen';
        }
      }
    }
  }"  
  
  search_hot_df = 
    reactiveVal(data.frame(`Mass` = rep("",10), `Ion` = "", `Mass Delta` = "", Term = "", Position = 0, check.names = FALSE))
  
  observe({
    search_df = 
      data.frame(Mass = unlist(str_split(input$search_hot_store, ",")), `Ion` = "", `Mass Delta` = "", Term = "", Position = 0, check.names = FALSE)
    search_hot_df(search_df)
  })
  
  observe({
    MASS_TOLERANCE = input$mass_tol_input
    OFF_BY_ONE = input$off_by_one_input
    
    if(is.null(input$search_hot)){
      return(data.frame(`Mass` = 0, `Ion` = "", `Mass Delta` = 0, Term = "", Position = 0, check.names = FALSE))
    }
    
    if(nrow(hot_to_r(input$search_hot)) == 0){
      return(data.frame(`Mass` = 0, `Ion` = "", `Mass Delta` = 0, Term = "", Position = 0, check.names = FALSE))
    }
    
    search_df = data.frame(Mass = hot_to_r(input$search_hot)$Mass, `Ion` = "", `Mass Delta` = 0, Term = "", Position = 0, check.names = FALSE)
    
    #Search
    ai = all_ions() %>% 
      mutate(mass = as.numeric(mass))
    
    #Clear past results
    search_df = search_df %>% 
      mutate(Ion = "") %>% 
      mutate(`Mass Delta` = "")
    
    rh = c()
    ch = c()
    
    
    for(i in 1:nrow(search_df)){
      result = search_df[i,] %>% 
        cross_join(ai) %>% 
        mutate(Mass = as.numeric(Mass)) %>% 
        mutate(`Mass Delta` = Mass - mass) %>% 
        mutate(`Mass Delta +1` = Mass - mass + 1) %>%
        mutate(`Mass Delta -1` = Mass - mass - 1) %>%
        mutate(`Ion` = ion_name) %>% 
        mutate(Term = term) %>% 
        mutate(Position = position) %>% 
        filter(abs(`Mass Delta`) < MASS_TOLERANCE | (OFF_BY_ONE & abs(abs(`Mass Delta`) - 1.0) < MASS_TOLERANCE)) %>% 
        rowwise() %>% 
        mutate(off = ifelse((OFF_BY_ONE & abs(abs(`Mass Delta`) - 1.0) < MASS_TOLERANCE), 
                            ifelse(abs(`Mass Delta +1`) < MASS_TOLERANCE, +1, -1), 0)) %>% 
        mutate(`Mass Delta` = ifelse(abs(`Mass Delta`) < MASS_TOLERANCE, `Mass Delta`,
                                     ifelse(abs(`Mass Delta +1`) < MASS_TOLERANCE, `Mass Delta +1`, `Mass Delta -1`))) %>% 
        mutate(`Mass Delta` = ifelse(input$mass_delta_input == "ppm", 
                                      formatC(`Mass Delta` / Mass * 1E6, digits = 2, format = "f"),
                                      formatC(`Mass Delta`, digits = ROUND_TO, format = "f")  )) %>% 
        mutate(`Mass Delta` = ifelse(off == 0, `Mass Delta`, paste0("(",`Mass Delta`, ")"))) %>% 
        ungroup() 
        
      if(nrow(result) > 0){
        if(input$table_view_input == "List"){
          rh = c(rh, result$list_row)
          ch = c(ch, result$list_col)
        }
        else {
          rh = c(rh, result$grid_row)
          ch = c(ch, result$grid_col)
        }
        
        search_df[i,] = result %>% 
          reframe(Mass, Ion = paste(Ion, collapse=", "), `Mass Delta` = paste(`Mass Delta`, collapse=", "), Term, Position)
      }
    }
    
    row_highlight(rh)
    col_highlight(ch)
    
    search_hot_df(search_df)
  })
  
  #Polarity button
  POLARITY = reactiveVal("+")
  observeEvent(input$polarity_input, {
    if(POLARITY() == "+")
      POLARITY("-")
    else if(POLARITY() == "-")
      POLARITY("+")
    updateActionButton(session, "polarity_input", label = POLARITY())
  })
  
  observeEvent(input$deuterium_exchange_input, {
    if(input$deuterium_exchange_input){
      ion_types_ref <<- ion_types_deut
      amino_acids_ref <<- amino_acids_deut
      term_ref <<- term_deut
    } 
    else{
      ion_types_ref <<- ion_types_norm
      amino_acids_ref <<- amino_acids_norm
      term_ref <<- term_norm
    }
  })
  

  
#Main code
  output$elemental_comp = renderText({
    sequence = sequence_input()
    Nterm = input$Nterm_input
    Cterm = input$Cterm_input
    deut = input$deuterium_exchange_input
    
    atoms = t(data.frame(sequence_to_atoms(sequence)))
    atoms = rbind(atoms, filter(term_ref, terminus == Nterm | terminus == Cterm) %>% select(C:D))
    atoms = atoms %>% summarise(across(C:D, ~ sum(.x, na.rm = TRUE)))
    
    label = ""
    for(a in 1:ncol(atoms)){
      if(atoms[,a] > 0)
        label = paste0(label, colnames(atoms)[a], "<sub>", atoms[,a], "</sub> ")
    }
    
    paste("Elemental Composition:", label)
  })
  
  output$aa_comp = renderText({
    sequence = sequence_input()
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
  
  output$manual_mods = renderText({
    mods = mods()
    paste("Manual modifications:", paste(mods$mod_mass, mods$mod_position, sep = " at #", collapse = ", "))
  })
   

  output$ion_hot = renderRHandsontable({
    req(row_highlight)
    req(col_highlight)
    
    ai = all_ions() %>% 
      mutate(mass = formatC(mass, digits = ROUND_TO, format = "f"))
    
    #Format for list or table
    if(input$table_view_input == "List"){
      ai = ai %>% 
        select(c("mass", "ion_name"))
    }
    else{
      ai = ai %>% 
        select(c("ion_type_charge", "mass", "length")) %>% 
        pivot_wider(names_from = length, values_from = mass) %>% 
        arrange(ion_type_charge)
    }
    
    rhandsontable(ai, col_highlight = col_highlight(), row_highlight = row_highlight(), 
                  rowHeaders = FALSE, readOnly = TRUE) %>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
      hot_cols(colWidths = 150) %>% 
      hot_cols(renderer = cell_renderer)
  })
  
  output$search_hot = renderRHandsontable({
    df = search_hot_df() %>% 
      select(Mass, Ion, 'Mass Delta')
    
    updateTextInput(session, "search_hot_store", value = paste0(df$Mass, collapse = ","))
    
    rhandsontable(df, rowHeaders = FALSE)%>%
      hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>% 
      hot_cols(renderer = htmlwidgets::JS("Handsontable.renderers.HtmlRenderer")) %>% 
      hot_col(2, readOnly = TRUE ) %>% 
      hot_col(3, readOnly = TRUE )
  })
  
  output$sequence_picture = renderUI({
    #  ┌ ┐
    df = search_hot_df()
    
    sequence = sequence_input()
    sequence = unlist(strsplit(sequence, ""))
    
    result = ""
    
    display_interval = 20
    
    for(a in 0:(length(sequence)/display_interval)){
      start = a*display_interval + 1
      end = a*display_interval + display_interval
      seq = sequence[start:end]
      
      seq_up = c("")
      seq_mid = c("")
      seq_down = c("")
      
      for(i in 1:length(seq)){
        seq_up[i*2] = " "
        if(!is.na(seq[i]))
          seq_mid[i*2] = seq[i]
        else
          seq_mid[i*2] = " "
        seq_down[i*2] = " "
        
        seq_up[i*2 + 1] = " "
        seq_mid[i*2 + 1] = " "
        seq_down[i*2 + 1] = " "
        
        search = filter(df, Position == i + start - 1)
        if(nrow(search) > 0){
          if("N" %in% search$Term){
            seq_up[i*2 + 1] = "┐"
          }
          if("C" %in% search$Term){
            seq_down[i*2 + 1] = "└"
          }
        }
      }
      
      #Add space at end
      #seq_up[i*2 + 1] = " "
      seq_mid[i*2 + 1] = paste0("   ", end)
      #seq_down[i*2 + 1] = " "
      
      seq_up = paste0(seq_up, collapse = "")
      seq_mid = paste0(seq_mid, collapse = "")
      seq_down = paste0(seq_down, collapse = "")
      
      result = paste(result, "<div style=\"font-family:'Courier New'\">", seq_up, "<br/>", seq_mid, "<br/>", seq_down, "</div>", sep = "", collapse = "")
    }
    
    HTML(result)
    #"ᒣᒥᒧᒪ"
  })
  
  
  #stores setup - insert at the bottom  !!!IMPORTANT
  appid = "application_fragmentor"
  setupStorage(appId = appid,inputs = TRUE)
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



