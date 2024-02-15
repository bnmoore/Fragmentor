library(tidyverse)
library(data.table)
library(Spectra)
library(readxl)
library(rbenchmark)


sequence_input = "MAAGVAAWLPFARAAAIGWMPVANCPMPLAPADKNKRQDELIVLNVSGRRFQTWRTTLERYPDTLLGSTEKEFFFNEDTKEYFFDRDPEVFRCVLNFYRTGKLHYPRYECISAYDDELAFYGILPEIIGDCCYEEYKDRKRENAERLMDDNDSENNQESMPSLSFRQTMWRAFENPHTSTLALVFYYVTGFFIAVSVITNVVETVPCGTVPGSKELPCGERYSVAFFCLDTACVMIFTVEYLLRLFAAPSRYRFIRSVMSIIDVVAIMPYYIGLVMTNNEDVSGAFVTLRVFRVFRIFKFSRHSQGLRILGYTLKSCASELGFLLFSLTMAIIIFATVMFYAEKGSSASKFTSIPASFWYTIVTMTTLGYGDMVPKTIAGKIFGSICSLSGVLVIALPVPVIVSNFSRIYHQNQRADKRRAQKKARLARIRVAKTGSSNAYLHSKRNGLLNEALELTGTPEEEHMGKTTSLIESQHHHLLHCLEKTTGLSYLVDDPLLSVRTSTIKNHEFIDEQMFEQNCMESSMQNYPSTRSPSLSSHPGLTTTCCSRRSKKTTHLPNSNLPATRLRSMQELSTIHIQGSEQPSLTTSRSSLNLKADDGLRPNCKTSQITTAIISIPTPPALTPEGESRPPPASPGPNTNIPSIASNVVKVSAL"

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

sequence_to_atoms = function(sequenceString, deut_exchange){
  aar = amino_acids_ref
  if(deut_exchange){
    aar = aar %>% 
      mutate(D = exchangeable_deuteriums) %>% 
      mutate(H = H - exchangeable_deuteriums)
  }
  
  aar = as.data.table(aar, key = "Abbrev1")
  atom_table = data.table(Abbrev1 = unlist(strsplit(sequenceString, split = "", fixed = TRUE)))
  atom_table = aar[match(atom_table$Abbrev1, aar$Abbrev1)]
  atom_table = atom_table[,!1:6]
  atom_table = colSums(atom_table)
  atom_table
}
