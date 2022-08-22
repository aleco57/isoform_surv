##################################################
#### TPM Extraction 
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
IDs <- read.csv(file = "Data/knownGene_ucsc_enst", header=F, sep=" ")

BRCA_match <- merge(BRCA_match, dplyr::filter(IDs, V2 %in% BRCA_match$Normal_Isoform), by.x="Normal_Isoform", by.y="V2", all=T) %>%
  dplyr::rename(Normal_ENST = V1)
BRCA_match <- merge(BRCA_match, dplyr::filter(IDs, V2 %in% BRCA_match$Cancer_Isoform), by.x="Cancer_Isoform", by.y="V2", all=T) %>%
  dplyr::rename(Cancer_ENST = V1)                    

BRCA_match[!complete.cases(BRCA_match),]
#Need to go back to try and match these missing ones, Sergio may be able to help


BRCA_match_cc <- na.omit(BRCA_match)

#Read in BRCA Survival data
BRCA_surv <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "survival/BRCA_survival.txt") %>% 
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()

##Read in tpsm data
load("/scratch/kvdd952/isoform_data/iso_ALL.RData")
#Subset to only BRCA survival data
sample <- BRCA_surv$sample %>% append("sample")
iso_BRCA <- iso_ALL[, colnames(iso_ALL) %in% sample]
#Remove after . so have the only unique identifier
iso_BRCA$sample <- gsub("\\..*","", iso_BRCA$sample)

#Filter to the specific isoform expression we want
Normal <- dplyr::filter(iso_BRCA, sample %in% BRCA_match_cc$Normal_ENST) %>% t() %>% janitor::row_to_names(1) 
Cancer <- dplyr::filter(iso_BRCA, sample %in% BRCA_match_cc$Cancer_ENST) %>% t() %>% janitor::row_to_names(1) 

#CHECK FOR MISSING!
## Write function to identify mismatch

#The value below is missing from  Healthy isoform expression:
# "ENST00000644486"

#Remove missing
BRCA_match_cc <- filter(BRCA_match_cc, Normal_ENST != "ENST00000644486")
#Repeat with cc
Normal <- dplyr::filter(iso_BRCA, sample %in% BRCA_match_cc$Normal_ENST) %>% t() %>% janitor::row_to_names(1) 
Cancer <- dplyr::filter(iso_BRCA, sample %in% BRCA_match_cc$Cancer_ENST) %>% t() %>% janitor::row_to_names(1) 

#Change our expression matrices to numeric
class(Normal) <- "numeric"
class(Cancer) <- "numeric"

#Match
Cancer <- Cancer[,match(BRCA_match_cc$Cancer_ENST, colnames(Cancer))]
Normal <- Normal[,match(BRCA_match_cc$Normal_ENST, colnames(Normal))]

#Save clean
save(BRCA_match, BRCA_match_cc, Cancer, Normal, file="Results/isoexp_clean")
