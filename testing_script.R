library(tidyverse)
library(data.table)
library(Spectra)
#library(rawrr)
library(readxl)


sequence_input = "RGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALG"
#file_path = "C:/Users/moore/Desktop/Test.raw"


#Load References
source("global.R")



atom_table_wide = atoms_ref %>% 
  select(Abbrev1, MonoisotopicMass) %>% 
  pivot_wider(names_from = "Abbrev1", values_from = "MonoisotopicMass")



# rawrr::installRawFileReaderDLLs()
# rawrr::installRawrrExe()
# rawrr::readSpectrum(rawfile = file_path, scan = 1)



#Parse sequence
Ntermbuild = substring(sequence_input, 1, 2:(nchar(sequence_input))-1)
Ctermbuild = substring(sequence_input, 2:(nchar(sequence_input)), nchar(sequence_input))


#Generate ions
all_ions = rbind(
  expand.grid(sequence = Ntermbuild, 
              ion_type = c("a", "b", "c"), 
              charge = seq(1, 3, 1), 
              terminus = "H", 
              stringsAsFactors = FALSE),
  expand.grid(sequence = Ctermbuild, 
              ion_type = c("x", "y", "z"), 
              charge = seq(1, 3, 1), 
              terminus = "OH", 
              stringsAsFactors = FALSE)
)

a0 = data.frame(t(sapply(unique(all_ions$sequence), sequence_to_atoms)), check.names = FALSE) %>% 
  mutate_all(as.numeric)
a0$sequence = rownames(a0)
all_ions0 = all_ions %>% 
  left_join(a0, by = "sequence")

all_ions1 = all_ions %>% 
  left_join(ion_types_ref, by = "ion_type")

all_ions2 = all_ions %>% 
  left_join(Nterm_ref, by = "terminus")

all_ions3 = all_ions %>%
  left_join(Cterm_ref, by = "terminus")

all_ions4 = all_ions %>%
  mutate(`H+` = charge)

all_ions98 = bind_rows(all_ions0, all_ions1, all_ions2, all_ions3, all_ions4)

all_ions = all_ions98 %>%
  group_by(sequence, ion_type, charge, terminus) %>% 
  summarise(across(C:`e-`, ~ sum(.x, na.rm = TRUE)))


all_ions = all_ions %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(ion_name = paste0(ion_type, str_length(sequence), " ", charge, "+")) %>% 
  mutate(ion_type_charge = paste0(ion_type, " ", charge, "+")) %>% 
  mutate(length = str_length(sequence))

#Mass calculation
mass = all_ions %>% 
  select(C:`e-`)
mass = as.matrix(mass) %*% diag(atom_table_wide)
mass = rowSums(mass)
all_ions = cbind(all_ions, mass)
all_ions = all_ions %>% 
  #mutate(mass = sum(c_across(C:`e-`) %*% atom_table_wide)) %>% 
  mutate(mass = mass / charge) %>% 
  arrange(mass) %>% 
  mutate(mass = format(mass, nsmall = 5))



sequenceString= "RGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALGRGYALG"
atom_table = data.frame(Abbrev1 = as.factor(unlist(strsplit(sequenceString, split = "")))) 
amino_acids_ref$Abbrev1 = as.factor(amino_acids_ref$Abbrev1)
atom_table = atom_table %>% 
  left_join(amino_acids_ref, by = "Abbrev1") %>% 
  select(C:`e-`)
atom_table = colSums(atom_table)
atom_table


