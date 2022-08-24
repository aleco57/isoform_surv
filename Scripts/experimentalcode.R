library(UCSCXenaTools)
library(biomaRt)
library(Homo.sapiens)

#Can use this a search tool:
#UCSCXenaTools::XenaShiny()

iso_ALL <- XenaGenerate() %>% 
   XenaFilter(filterDatasets = "TcgaTargetGtex_rsem_isoform_tpm") %>%
   XenaQuery() %>%
   XenaDownload() %>%
   XenaPrepare()

iso_pct <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "TcgaTargetGtex_rsem_isopct") %>%
  XenaQuery() %>%
  XenaDownload()
  

l <- BRCA_surv[["sample"]]
l <- append(l, "sample")
iso_pct <- data.table::fread( file="/scratch/kvdd952/isoform_data/TcgaTargetGtex_rsem_isopct", sep = "\t", select=l)

save(iso_pct, file="/scratch/kvdd952/isoform_data/iso_pct_BRCA.RData")

#"TcgaTargetGtex_rsem_isopct"


#Filter the R data to just the BRCA patients we are interested in
BRCA_surv <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "survival/BRCA_survival.txt") %>% 
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()

load("/scratch/kvdd952/isoform_data/iso_ALL.RData")

#Subset to only BRCA survival data
sample <- BRCA_surv$sample %>% append("sample")
iso_BRCA <- iso_ALL[, colnames(iso_ALL) %in% sample]


#Read in USCS_ENST converter
IDs <- read.csv(file = "Data/UCSC_ENST.txt", sep="\t", header=F)

BRCA_normal <- c("uc003gbj", "uc004aso", "uc003ncb")
BRCA_cancer <- c("uc003gbi", "uc004asp", "uc011djb")

BRCA_normalENST<- dplyr::filter(IDs, grepl(paste(BRCA_normal, collapse="|"), V1)) [["V2"]]
BRCA_cancerENST <- dplyr::filter(IDs, grepl(paste(BRCA_cancer, collapse="|"), V1)) [["V2"]]

BRCA_normalVar <- dplyr::filter(iso_BRCA, grepl(paste(BRCA_normalENST, collapse="|"), sample)) 
BRCA_cancerVar <- dplyr::filter(iso_BRCA, grepl(paste(BRCA_cancerENST, collapse="|"), sample))

C <- t(BRCA_cancerVar) %>% janitor::row_to_names(1)
C <- as.numeric(C)
N <- t(BRCA_normalVar) %>% janitor::row_to_names(1) %>% as.data.frame()
N <- sapply(N, as.numeric)

CLPX1 <- C[,1] - N[,1]
BICD2 <- C[,2] - N[,2]
CAP2 <-  C[,3] - N[,3]

#Subset so only pan cancer switches for BRCA
#identify the isoforms on ensemble
HR <- readxl::read_xlsx("~/Documents/SupplementaryTab_IsoSwitch/Table_S5_HR.xlsx", skip=1)
#Identify transcript variant through Biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

dplyr::filter(HR, BRCA >5) [["gene_name"]]
dplyr::filter(iso_BRCA, grepl("ENSG00000136717"))

cli %>% dplyr::filter(grepl("ENSG00000242268", Ensembl_ID))


library(Homo.sapiens)
res1 <- select(Homo.sapiens,
               keys(Homo.sapiens, keytype="TXID"),
               columns=c("GENEID","TXNAME","TXCHROM", "ENSEMBLTRANS"), keytype="TXID")

dplyr::filter(res1, grepl(paste(BRCA_match$Normal_Isoform, collapse="|"), TXNAME))


#Function to identify mismatch:
for (i in c("Cancer", "Normal")) {
  if (BRCA_match_cc[paste(i, "_ENST", sep="")] == colnames(i)) {
    print("MATCH")
  } else {
    setdiff(BRCA_match_cc[paste(i, "_ENST", sep="")], colnames(i))
  } 
  
}
nrow(BRCA_match_cc)
ncol(Normal_exp)
ncol(Cancer_exp)

setdiff(BRCA_match_cc$Normal_ENST, colnames(Normal_exp))
