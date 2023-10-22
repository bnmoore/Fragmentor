library(stringr)
library(tidyverse)
library(readxl)
library(data.table)
library(BiocManager)
#library(rawrr)

#Load References
ion_types_ref = read_excel("resources.xlsx", sheet = "ion_type_list")
atoms_ref = read_excel("resources.xlsx", sheet = "atoms") #Atomic masses from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
amino_acids_ref = read_excel("resources.xlsx", sheet = "amino_acids")
losses_ref = read_excel("resources.xlsx", sheet = "losses")
term_ref = read_excel("resources.xlsx", sheet = "term")


#Atom table
atom_table_wide_mono = atoms_ref %>% 
  select(Abbrev1, MonoisotopicMass) %>% 
  pivot_wider(names_from = "Abbrev1", values_from = "MonoisotopicMass")

atom_table_wide_avg = atoms_ref %>% 
  select(Abbrev1, AverageMass) %>% 
  pivot_wider(names_from = "Abbrev1", values_from = "AverageMass")


#Sequence to atoms
sequence_to_atoms = function(sequenceString, deut_exchange){
  aar = amino_acids_ref
  if(deut_exchange){
    aar = aar %>% 
      mutate(D = exchangeable_deuteriums) %>% 
      mutate(H = H - exchangeable_deuteriums)
  }
  
  atom_table = data.frame(Abbrev1 = unlist(strsplit(sequenceString, split = ""))) 
  atom_table = atom_table %>% 
    left_join(aar, by = "Abbrev1") %>% 
    select(C:D)
  atom_table = colSums(atom_table)
  atom_table
}




#Atoms to mass
atoms_to_mass = function(atom_table, mono){
  if(mono == "Monoisotopic"){
    for(i in 1:length(atoms_ref$Abbrev1)){
      atom_table[1,i] = atom_table[1,i] * atoms_ref$MonoisotopicMass[i]
    }
  }
  else{
    for(i in 1:length(atoms_ref$Abbrev1)){
      atom_table[1,i] = atom_table[1,i] * atoms_ref$AverageMass[i]
    }
  }
  sum(atom_table)
}

#Naming for losses
losses_ref = losses_ref %>% 
  rowwise() %>% 
  mutate(loss_mass = atoms_to_mass(across(C:D), TRUE)) %>% 
  mutate(loss_display = if_else(str_detect(loss, "sidechain"), 
                                   paste0("sidechain"),
                                   paste0(loss))) %>% 
  mutate(loss = if_else(str_detect(loss, "sidechain"), 
                                paste0(round(loss_mass), Abbrev1),
                                paste0(loss)))








# #Calculate amino acid masses
# amino_acids_ref = amino_acids_ref %>% 
#   rowwise() %>% 
#   mutate(MonoisotopicMass = atoms_to_mass(across(C:`e-`), TRUE)) %>% 
#   mutate(AverageMass = atoms_to_mass(across(C:`e-`), FALSE))
# 
# #Calculate ion type masses
# ion_types_ref = ion_types_ref %>% 
#   rowwise() %>% 
#   mutate(MonoisotopicMass = atoms_to_mass(across(C:`e-`), TRUE)) %>% 
#   mutate(AverageMass = atoms_to_mass(across(C:`e-`), FALSE))
#
# #Atoms to mass
# atoms_to_mass = function(atom_table){
#   mass_table = atom_table %>% 
#     pivot_longer(cols = everything(), names_to = "Abbrev1", values_to = "number") 
#   mass_table = mass_table %>%  
#     left_join(atoms_ref, by = "Abbrev1") 
#   mass_table = mass_table %>% 
#     mutate(MonoisotopicMass = number * MonoisotopicMass) %>% 
#     mutate(AverageMass = number * AverageMass) %>% 
#     summarise(across(MonoisotopicMass:AverageMass, sum))
#   
#   mass_table
# }
