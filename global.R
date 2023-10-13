library(stringr)
library(tidyverse)
library(readxl)
library(data.table)

ionTypeList = read_excel("resources.xlsx", sheet = "ion_type_list")

#Atomic masses from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
atoms = read_excel("resources.xlsx", sheet = "atoms")

aminoAcids = read_excel("resources.xlsx", sheet = "amino_acids")
aminoAcids = data.table(aminoAcids)

losses = read_excel("resources.xlsx", sheet = "losses")

