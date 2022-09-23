##################################################
#### TPM Extraction 
##################################################

############
# BRCA as proof of concept
############

library(UCSCXenaTools)
library(tidyverse)
library(Homo.sapiens)

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
ID2 <- read.csv(file = "Data/UCSC_ENST.txt", header=F, sep="\t")
ID2$V1 <- gsub("\\..*","", ID2$V1)

BRCA_match <- merge(BRCA_match, dplyr::filter(ID2, V1 %in% BRCA_match$Cancer_Isoform), by.x="Cancer_Isoform", by.y="V1", all=T) %>%
  dplyr::rename(Cancer_ENST = V2)    

#Fill in missing isoforms from online search

#Label Genes
k<-keys(Homo.sapiens,keytype="TXNAME")
map <- select(Homo.sapiens,keys=k,columns=c("TXNAME","SYMBOL"),keytype="TXNAME")
map$TXNAME <- gsub("\\..*","", map$TXNAME)
BRCA_match <- merge(BRCA_match, map, by.x = "Cancer_Isoform", by.y = "TXNAME", all.x=T)

#There are four genes missing, three of these can be imputed by the normal isoform id
map %>% filter(TXNAME %in% BRCA_match[is.na(BRCA_match$SYMBOL),"Normal_Isoform"])
BRCA_match$SYMBOL[BRCA_match$Normal_Isoform == "uc010ogc"] <- "NKAIN1"
BRCA_match$SYMBOL[BRCA_match$Normal_Isoform == "uc003vej"] <- "BCAP29"
BRCA_match$SYMBOL[BRCA_match$Normal_Isoform == "uc001sws"] <- "TMEM19"

#The last has been imputed by manual search on Google
BRCA_match$SYMBOL[BRCA_match$Normal_Isoform == "uc010thn"] <- "DACH1"

#Add a the gene label for further analyses
BRCA_match$CancerIso <- paste(BRCA_match$SYMBOL,"CI", sep ="_")

#Drop the missing cancer isoforms
BRCA_match_cc <- dplyr::filter(BRCA_match, !is.na(Cancer_ENST))

#Read in BRCA Survival data
BRCA_surv <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "survival/BRCA_survival.txt") %>% 
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()

#Filter to IDs ending in "01" meaning all tumour samples
BRCA_surv <- BRCA_surv[endsWith(BRCA_surv$sample, '-01'),]

#Read in gene expression data
g <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "TcgaTargetGtex_RSEM_Hugo_norm_count") %>%
  XenaQuery()

Gene_exp <- fetch_dense_values(host = g$hosts,
                               dataset = g$datasets,
                               identifiers = BRCA_match_cc[["SYMBOL"]],
                               samples = BRCA_surv[["sample"]]) %>% t()


load("/scratch/kvdd952/isoform_data/iso_pct_BRCA.RData") #issues reading whole thing in, this is restricted to BRCA survival data patients
#Subset to only BRCA survival data
sample <- BRCA_surv$sample %>% append("sample")
iso_pct <- iso_pct %>% dplyr::select(matches(sample))
#Remove after . so have the only unique identifier
iso_pct$sample <- gsub("\\..*","", iso_pct$sample)
iso_pct <- iso_pct %>% relocate(sample)


#Repeat with iso%
Cancer_pct <- dplyr::filter(iso_pct, sample %in% BRCA_match_cc$Cancer_ENST) %>% t() %>% janitor::row_to_names(1)
class(Cancer_pct) <- "numeric"

#Restrict analyses to only those with isoform expression available, (four isoforms missing)
BRCA_match_cc <- BRCA_match_cc %>% dplyr::filter(Cancer_ENST %in% colnames(Cancer_pct))

#Which way does this need to be matched to change columns in Cancer_pct?
Cancer_pct<-Cancer_pct[,match(BRCA_match_cc$Cancer_ENST, colnames(Cancer_pct))]
identical(colnames(Cancer_pct),BRCA_match_cc$Cancer_ENST) %>% stopifnot()
colnames(Cancer_pct) <- BRCA_match_cc$CancerIso

### Extracting Gene Expression
##############################


#Now need to extract gene information from Xena
g <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "TcgaTargetGtex_RSEM_Hugo_norm_count") %>%
  XenaQuery()

Gene_exp <- fetch_dense_values(host = g@hosts,
                               dataset = g@datasets,
                               identifiers = BRCA_match[["SYMBOL"]],
                               samples = row.names(BRCA_surv)) %>% t()

#Read in Phenotype
BRCA_phen <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") %>% 
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()

#Filter to only patients available
BRCA_phen <- dplyr::filter(BRCA_phen, sampleID %in% row.names(Cancer_pct)) %>% column_to_rownames(var="sampleID")

#Mutate and split stage into five categories
BRCA_phen <- mutate(BRCA_phen, stage = case_when(
  pathologic_stage == "Stage I" | pathologic_stage == "Stage IA" | pathologic_stage == "Stage IB" ~ 1,
  pathologic_stage == "Stage II" | pathologic_stage == "Stage IIA" | pathologic_stage == "Stage IIB" ~ 2,
  pathologic_stage == "Stage III" | pathologic_stage == "Stage IIIA" | pathologic_stage == "Stage IIIB" | pathologic_stage == "Stage IIIC" ~ 3,
  pathologic_stage == "Stage IV" ~ 4,
  pathologic_stage == "Stage X" ~ 5)) 
BRCA_phen$stage <- as.factor(BRCA_phen$stage)
BRCA_phen$gender <- ifelse(BRCA_phen$gender=="FEMALE", 0, 1)

#Merge
clean_data_all <- merge(Cancer_pct, BRCA_surv[,c(1,3,4)], by.x=0, by.y="sample") %>% column_to_rownames(var="Row.names") %>% relocate(OS,OS.time)
clean_data_all <- merge(clean_data_all, Gene_exp, by=0) %>% column_to_rownames(var="Row.names")
clean_data_all <- merge(clean_data_all, BRCA_phen[, c("stage", "gender", "age_at_initial_pathologic_diagnosis", "days_to_birth")], by=0) %>% column_to_rownames(var="Row.names")


#Remove those if don't have data on both gene and isoform expression
clean_data_all <- subset(clean_data_all, select=-c(CENATAC_CI, EFL1_CI, DMAC2_CI, 
                                 MAP3K20_CI, GRAMD2A_CI, 
                                 CLPTM1L, TSPAN17, ERBB2, FAM228B))

save(clean_data_all, file = "Results/clean_data_all.RData")

surv_data <- 
surv_data[is.na(surv_data$OS.time),]
missing <- dplyr::filter(BRCA_surv, !sample %in% row.names(surv_data))
missing <- missing[endsWith(missing$sample, '-01'),]



merge(surv_data, BRCA_phen[, c("stage", "gender", "days_to_birth")], by=0) 




#Previous Code
#IDs <- read.csv(file = "Data/knownGene_ucsc_enst", header=F, sep=" ")
#There may be missing matches in this dataframe below

#BRCA_match <- merge(BRCA_match, dplyr::filter(IDs, V2 %in% BRCA_match$Normal_Isoform), by.x="Normal_Isoform", by.y="V2", all=T) %>%
dplyr::rename(Normal_ENST = V1)
#BRCA_match <- merge(BRCA_match, dplyr::filter(IDs, V2 %in% BRCA_match$Cancer_Isoform), by.x="Cancer_Isoform", by.y="V2", all=T) %>%
dplyr::rename(Cancer_ENST = V1)                    

BRCA_match[!complete.cases(BRCA_match),]
#Need to go back to try and match these missing ones

##Read in %iso data
#load("/scratch/kvdd952/isoform_data/iso_ALL.RData") #This is TPM data which have ignored for now
