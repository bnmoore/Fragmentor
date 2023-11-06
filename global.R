library(stringr)
library(tidyverse)
library(readxl)
library(data.table)
#library(BiocManager)
#library(rawrr)

ROUND_TO = 5

#Load References
ion_types_norm = read_excel("resources.xlsx", sheet = "ion_type_norm")
ion_types_deut = read_excel("resources.xlsx", sheet = "ion_type_deut")
atoms_ref = read_excel("resources.xlsx", sheet = "atoms") #Atomic masses from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
amino_acids_norm = read_excel("resources.xlsx", sheet = "amino_acids_norm")
amino_acids_deut = read_excel("resources.xlsx", sheet = "amino_acids_deut")
losses_ref = read_excel("resources.xlsx", sheet = "losses")
term_norm = read_excel("resources.xlsx", sheet = "term_norm")
term_deut = read_excel("resources.xlsx", sheet = "term_deut")

#Global tables
ion_types_ref = ion_types_norm
amino_acids_ref = amino_acids_norm
term_ref = term_norm

#Atom table
atom_table_wide_mono = atoms_ref %>% 
  select(Abbrev1, MonoisotopicMass) %>% 
  pivot_wider(names_from = "Abbrev1", values_from = "MonoisotopicMass")

atom_table_wide_avg = atoms_ref %>% 
  select(Abbrev1, AverageMass) %>% 
  pivot_wider(names_from = "Abbrev1", values_from = "AverageMass")

#Amino acids norm
amino_acids_norm = amino_acids_norm %>% 
  arrange(Abbrev1)
amino_acids_norm = as.data.table(amino_acids_norm, key = "Abbrev1")

#Amino acids deut
amino_acids_deut = amino_acids_deut %>% 
  arrange(Abbrev1)
amino_acids_deut = as.data.table(amino_acids_deut, key = "Abbrev1")

#Sequence to atoms
sequence_to_atoms = function(sequenceString){
  aar = amino_acids_ref
  atom_table = data.table(Abbrev1 = unlist(strsplit(sequenceString, split = "", fixed = TRUE)))
  atom_table = aar[match(atom_table$Abbrev1, aar$Abbrev1)]
  atom_table = atom_table[,!1:5]
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
  mutate(loss_display = if_else(sidechain == 1, paste0(round(loss_mass), Abbrev1), loss))







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
