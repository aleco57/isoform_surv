##################################################
#### Pipeline for isoform survival analysis ######
##################################################

############
# BRCA as proof of concept
############

library(UCSCXenaTools)
library(tidyverse)


# (1) Download and extract relevant isoform switches
Clemente <- read.csv(file = "Data/IsoSwitch_meta/Clemente/BRCA_cancer", sep="_")
Kahraman <- read.csv(file = "Data/IsoSwitch_meta/Kahraman/BRCA_cancer", sep="_")
Vitting <- read.csv(file = "Data/IsoSwitch_meta/Vitting/BRCA_cancer", sep="_")

BRCA_match <- inner_join(Clemente, Kahraman) %>% 
  rbind(inner_join(Clemente, Vitting)) %>% 
  rbind(inner_join(Vitting, Kahraman)) %>% unique()

#Identify those switches found in all three studies
BRCA_match$tripmatch <- ifelse(do.call(paste0, BRCA_match) %in% do.call(paste0, inner_join(Clemente, inner_join(Vitting, Kahraman))),1,0)

#Extract UCSC to ENSEMBL ID Converter
IDs <- read.csv(file = "Data/UCSC_ENST.txt", sep="\t", header=F)
