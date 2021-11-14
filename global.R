library(stringr)
library(tidyverse)
library(OrgMassSpecR)
library(data.table)

ionTypeList = list("M", "M*", "a", "b", "b+H2O", "c", "x", "y", "z", "d", "w")

#Atomic masses from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
atoms = read.csv("atoms.csv", stringsAsFactors = FALSE)

aminoAcids = read.csv("amino_acids.csv", stringsAsFactors = FALSE)
aminoAcids$MonoisotopicMass = sapply(aminoAcids$Composition, function(x){ MonoisotopicMass(formula = ListFormula(x))})
aminoAcids$MolecularWeight = sapply(aminoAcids$Composition, function(x){ MolecularWeight(formula = ListFormula(x))})
aminoAcids = data.table(aminoAcids, key = "Abbrev1")

losses = read.csv("losses.csv", stringsAsFactors = FALSE)
losses[is.na(losses)] = ""
losses$MonoisotopicMass = sapply(losses$Gain, function(x){ MonoisotopicMass(formula = ListFormula(x))})
losses$MolecularWeight = sapply(losses$Gain, function(x){ MolecularWeight(formula = ListFormula(x))})
losses$MonoisotopicMass = losses$MonoisotopicMass - sapply(losses$Loss, function(x){ MonoisotopicMass(formula = ListFormula(x))})
losses$MolecularWeight = losses$MolecularWeight - sapply(losses$Loss, function(x){ MolecularWeight(formula = ListFormula(x))})