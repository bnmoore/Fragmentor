library(tidyverse)
library(data.table)
library(Spectra)
#library(rawrr)
library(readxl)
library(fuzzyjoin)
library(rbenchmark)


sequence_input = "MAAGVAAWLPFARAAAIGWMPVANCPMPLAPADKNKRQDELIVLNVSGRRFQTWRTTLERYPDTLLGSTEKEFFFNEDTKEYFFDRDPEVFRCVLNFYRTGKLHYPRYECISAYDDELAFYGILPEIIGDCCYEEYKDRKRENAERLMDDNDSENNQESMPSLSFRQTMWRAFENPHTSTLALVFYYVTGFFIAVSVITNVVETVPCGTVPGSKELPCGERYSVAFFCLDTACVMIFTVEYLLRLFAAPSRYRFIRSVMSIIDVVAIMPYYIGLVMTNNEDVSGAFVTLRVFRVFRIFKFSRHSQGLRILGYTLKSCASELGFLLFSLTMAIIIFATVMFYAEKGSSASKFTSIPASFWYTIVTMTTLGYGDMVPKTIAGKIFGSICSLSGVLVIALPVPVIVSNFSRIYHQNQRADKRRAQKKARLARIRVAKTGSSNAYLHSKRNGLLNEALELTGTPEEEHMGKTTSLIESQHHHLLHCLEKTTGLSYLVDDPLLSVRTSTIKNHEFIDEQMFEQNCMESSMQNYPSTRSPSLSSHPGLTTTCCSRRSKKTTHLPNSNLPATRLRSMQELSTIHIQGSEQPSLTTSRSSLNLKADDGLRPNCKTSQITTAIISIPTPPALTPEGESRPPPASPGPNTNIPSIASNVVKVSAL"
#file_path = "C:/Users/moore/Desktop/Test.raw"

#Sequence to atoms
benchmark(s2a = {
  aar = as.data.table(amino_acids_ref, key = "Abbrev1")
  atom_table = data.table(Abbrev1 = unlist(strsplit(sequence_input, split = "")), key = "Abbrev1")
  #atom_table = aar[atom_table, on = 'Abbrev1']
  atom_table = aar[match(atom_table$Abbrev1, aar$Abbrev1)]
  atom_table = atom_table[,!1:6]
  #atom_table = atom_table[, lapply(.SD, sum)]
  #atom_table = atom_table[, colSums(.SD)]
  atom_table = colSums(atom_table)
  atom_table
  }
)

#Load References
source("global.R")


# rawrr::installRawFileReaderDLLs()
# rawrr::installRawrrExe()
# rawrr::readSpectrum(rawfile = file_path, scan = 1)



#Parse sequence
Mbuild = sequence_input
Ntermbuild = substring(sequence_input, 1, 2:(nchar(sequence_input))-1)
Ctermbuild = substring(sequence_input, 2:(nchar(sequence_input)), nchar(sequence_input))


#Generate ions
all_ions = rbind(
  expand.grid(sequence = Mbuild, 
              ion_type = "M", 
              charge = seq(1, 3, 1), 
              terminus = "M", 
              loss = c("", "NH3", "H2O"),
              stringsAsFactors = FALSE),
  expand.grid(sequence = Ntermbuild, 
              ion_type = "d", 
              charge = seq(1, 3, 1), 
              terminus = "H", 
              loss = c("", "NH3", "H2O"),
              stringsAsFactors = FALSE),
  expand.grid(sequence = Ctermbuild, 
              ion_type = "y", 
              charge = seq(1, 3, 1), 
              terminus = "OH", 
              loss = c("", "NH3", "H2O"),
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
  mutate(terminus = "H")%>% 
  left_join(term_ref, by = "terminus") 

all_ions4 = all_ions %>% 
  filter(terminus == "M") %>% 
  mutate(terminus = "OH") %>% 
  left_join(term_ref, by = "terminus") 

all_ions5 = all_ions %>%
  mutate(`H+` = charge)

all_ions6 = all_ions %>% 
  filter(terminus != "M") %>% 
  left_join(term_ref, by = "terminus") %>% 
  mutate(last_aa = if_else(term == "N", substring(sequence, nchar(sequence), nchar(sequence)), "")) %>% 
  mutate(last_aa = if_else(term == "C", substring(sequence, 1, 1), last_aa)) %>% 
  select(-(C:D)) %>% 
  filter(ion_type == "d") %>% 
  left_join(filter(losses_ref, loss == "partial_sidechain"), 
            by = c("last_aa" = "Abbrev1"), relationship = "many-to-many")


all_ions7 = all_ions %>% 
  left_join(losses_ref, by = "loss", relationship = "many-to-many") %>% 
  group_by(sequence, ion_type, charge, loss) %>% 
  filter(row_number() == 1)



all_ions98 = bind_rows(all_ions0, all_ions1, all_ions2, all_ions3, all_ions4, all_ions5, all_ions6, all_ions7)

all_ions = all_ions98 %>%
  ungroup() %>% 
  group_by(sequence, ion_type, charge, loss) %>% 
  summarise(across(C:D, ~ sum(.x, na.rm = TRUE)))

polarity = "+"

all_ions = all_ions %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(ion_name = paste0(ion_type, "<sub>", str_length(sequence), "</sub> ", charge, polarity)) %>% 
  mutate(ion_type_charge = paste0(ion_type, " ", charge, polarity)) %>% 
  mutate(length = str_length(sequence))

#Mass calculation
mass = all_ions %>% 
  select(C:D)
mass = as.matrix(mass) %*% diag(atom_table_wide)
mass = rowSums(mass)
all_ions = cbind(all_ions, mass)
all_ions = all_ions %>% 
  mutate(mass = mass / charge) %>% 
  arrange(mass) %>% 
  mutate(mass = format(mass, nsmall = 5))




