############ Analysis Pipeline ############
################################################


library(UCSCXenaTools)
library(tidyverse)
library(Homo.sapiens)
library(corrr)
library(survival)
library(xgboost)
library(caret)
library(SHAPforxgboost)
library(pROC)
library(cluster)
library(factoextra)
library(janitor)
library(klaR)
library(nnet)
library(survminer)



#### (1) Extract TCGA data through Xena

#Read in switches from Github

dir <- "https://raw.githubusercontent.com/KarakulakTulay/Isoform_Comparison/master/Overlaps/IsoSwitch"

Clemente <- read.csv(file = file.path(dir, "Clemente/BRCA_cancer"), sep="_")
Kahraman  <- read.csv(file = file.path(dir, "Kahraman/BRCA_cancer"), sep="_")             
Vitting <- read.csv(file = file.path(dir, "Vitting/BRCA_cancer"), sep="_")

#Create a dataframe of those switches found in a minimum of two studies
BRCA_match <- inner_join(Clemente, Kahraman) %>% 
  rbind(inner_join(Clemente, Vitting)) %>% 
  rbind(inner_join(Vitting, Kahraman)) %>% unique()

#Around 5000 switches found in all the studies:
switches <- rbind(Clemente, Kahraman, Vitting) %>% unique()


### Convert UCSC naming to ENSEMBL
ID <- read.csv(file = "/Users/am17168/Downloads/UCSC_ENST.txt",
              header=F, sep="\t")
ID$V1 <- gsub("\\..*","", ID$V1)
BRCA_match <- merge(BRCA_match, dplyr::filter(ID, V1 %in% BRCA_match$Cancer_Isoform), by.x="Cancer_Isoform", by.y="V1", all=T) %>%
  dplyr::rename(Cancer_ENST = V2)


#Label isoforms with gene names
k<-keys(Homo.sapiens,keytype="TXNAME")
map <- AnnotationDbi::select(Homo.sapiens,keys=k,columns=c("TXNAME","SYMBOL"),keytype="TXNAME")
map$TXNAME <- gsub("\\..*","", map$TXNAME)
BRCA_match <- merge(BRCA_match, map, by.x = "Cancer_Isoform", by.y = "TXNAME", all.x=T)

#Add a the gene label for further analyses
BRCA_match$CancerIso <- paste(BRCA_match$SYMBOL,"CI", sep ="_")

#Drop the missing cancer isoforms
BRCA_match_cc <- dplyr::filter(BRCA_match, !is.na(Cancer_ENST))




########## Read in TCGA data from Xena

### Read in phenotype data
BRCA_phen <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") %>% 
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()

BRCA_surv <- XenaGenerate() %>% 
  XenaFilter(filterDatasets = "survival/BRCA_survival.txt") %>% 
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()

#Only 01 which is cancer tissue
BRCA_phen <- BRCA_phen[(grepl("-01", BRCA_phen$sampleID)),] 
BRCA_surv <- BRCA_surv[(grepl("-01", BRCA_surv$sample)),]


#Extract TPM data

all_identifiers <- fetch_dataset_identifiers("https://toil.xenahubs.net", 
                                             "TcgaTargetGtex_expected_count")
iso_pct <- fetch_dense_values("https://toil.xenahubs.net", dataset = "TcgaTargetGtex_rsem_isopct", samples = BRCA_phen$sampleID,
                               identifiers = all_identifiers[pmatch(BRCA_match_cc[["Cancer_ENST"]], all_identifiers)],
                               check = T)
iso_tpm <-  fetch_dense_values("https://toil.xenahubs.net", dataset = "TcgaTargetGtex_rsem_isoform_tpm", samples = BRCA_phen$sampleID,
                              identifiers = all_identifiers[pmatch(BRCA_match_cc[["Cancer_ENST"]], all_identifiers)],
                              check = T)
rownames(iso_pct) <- gsub("\\..*","", rownames(iso_pct))
rownames(iso_tpm) <- gsub("\\..*","", rownames(iso_tpm))

### "TcgaTargetGtex_expected_count", TcgaTargetGtex_rsem_isoform_tpm Use this to extract the expected count over to illustrate correlated variables

#Extract Gene expression data
Gene_exp <- fetch_dense_values(host = "https://toil.xenahubs.net",
                               dataset = "TcgaTargetGtex_RSEM_Hugo_norm_count",
                               identifiers = BRCA_match_cc[["SYMBOL"]],
                               samples = BRCA_phen$sampleID,
                               check=T)

#Add columns if gene/iso values present
BRCA_match_cc$Gene_exp.has <- ifelse(BRCA_match_cc$SYMBOL %in% rownames(Gene_exp), T, F)
BRCA_match_cc$Iso_exp.has <- ifelse(BRCA_match_cc$Cancer_ENST %in% rownames(iso_pct), T, F)

#Generate a stage variable from phenotype data
BRCA_phen <- mutate(BRCA_phen, stage = case_when(
  pathologic_stage == "Stage I" | pathologic_stage == "Stage IA" | pathologic_stage == "Stage IB" ~ 1,
  pathologic_stage == "Stage II" | pathologic_stage == "Stage IIA" | pathologic_stage == "Stage IIB" ~ 2,
  pathologic_stage == "Stage III" | pathologic_stage == "Stage IIIA" | pathologic_stage == "Stage IIIB" | pathologic_stage == "Stage IIIC" ~ 3,
  pathologic_stage == "Stage IV" ~ 4,
  pathologic_stage == "Stage X" ~ 5)) 

#Add a binary stage variable
BRCA_phen <- BRCA_phen %>% mutate(stage = na_if(stage, 5)) %>%
  mutate(stage_binary = case_when(stage == 1 | stage == 2 ~ 0,
                                  stage == 3 | stage == 4 ~ 1))


#Change age variable
BRCA_phen$age_days <- BRCA_phen$days_to_birth *-1

## Change ENST names to gene_CI (cancer isoform)
identical(filter(BRCA_match_cc, Iso_exp.has)[["Cancer_ENST"]], rownames(iso_pct)) %>% stopifnot()
rownames(iso_pct) <- filter(BRCA_match_cc, Iso_exp.has) [["CancerIso"]]

identical(filter(BRCA_match_cc, Iso_exp.has)[["Cancer_ENST"]], rownames(iso_tpm)) %>% stopifnot()
rownames(iso_tpm) <- filter(BRCA_match_cc, Iso_exp.has) [["CancerIso"]]

#Remove gene expression if we don't have isoform expression value
Gene_exp <- Gene_exp[(rownames(Gene_exp) %in% filter(BRCA_match_cc, Iso_exp.has) [["SYMBOL"]] %>% which()) ,]

## Clean Phenotype data and merge with survival + expression data
phen_exp <- BRCA_phen %>% dplyr::select(sampleID,
                                        PAM50Call_RNAseq, #Least NA in PAM50 variables
                                        #gender, #this variable should be extracted for other cancers
                                        stage,
                                        stage_binary,
                                        age_at_initial_pathologic_diagnosis,
                                        age_days) %>% 
  merge(BRCA_surv[,c(1,3,4)], by.x = "sampleID", by.y = "sample", all=T) %>% 
  merge(t(iso_pct), by.x = "sampleID", by.y = 0, all=T) %>%
  merge(t(Gene_exp), by.x = "sampleID", by.y = 0, all=T) %>%
  column_to_rownames(var = "sampleID")

# Change class of variables
phen_exp[, c("PAM50Call_RNAseq", "stage", "OS")] <- 
  sapply(phen_exp[, c("PAM50Call_RNAseq", "stage", "OS")] , factor)

phen_exp[, !colnames(phen_exp) %in% c("PAM50Call_RNAseq", "stage", "OS")] <- 
  sapply(phen_exp[, !colnames(phen_exp) %in% c("PAM50Call_RNAseq", "stage", "OS")] , as.numeric)

#Remove variables with a low dynamic range
range1 <- function(x){
  x <- na.omit(x)
  max(x) - min(x)
}
ten <- vector()
twenty <- vector()

for (col in colnames(phen_exp %>% dplyr::select(ends_with("_ci")))) {
  if (range1(phen_exp[col]) < 10) {
    ten <- c(ten, col)
  }
  if (range1(phen_exp[col]) < 20) {
    twenty <- c(ten, col)
  }
}

print(ten)
print(twenty)

phen_exp <- phen_exp[,-(match(sub("\\_.*", "", twenty), colnames(phen_exp)))]
phen_exp <- phen_exp[,-(match(twenty, colnames(phen_exp)))]

phen_exp <- janitor::clean_names(phen_exp)
colnames(phen_exp)[1] <- "pam50"

##### (2) EDA
# res.cor <- correlate(phen_exp %>% dplyr::select(-os, -os_time, -stage_binary, -stage, -age_days, -age_at_initial_pathologic_diagnosis, -pam50), method= "spearman")

#p <- res.cor %>%  
#  gather(-term, key = "colname", value = "cor") %>% 
#  filter(abs(cor) > 0.7) %>%
#  filter(!str_detect(term, "_ci"))

pct_cor <- vector()
tpm_cor <- vector()

for (i in 1:nrow(filter(BRCA_match_cc, Gene_exp.has, Iso_exp.has))) {
data <- filter(BRCA_match_cc, Gene_exp.has, Iso_exp.has)[i, c("SYMBOL", "CancerIso")]
pct_cor <- c(pct_cor, cor(iso_pct[data[[2]],], Gene_exp[data[[1]],], method = "spearman"))
tpm_cor <- c(tpm_cor, cor(iso_tpm[data[[2]],], Gene_exp[data[[1]],], method = "spearman"))
}

correlation_df <- tibble(Gene = filter(BRCA_match_cc, Gene_exp.has, Iso_exp.has)[["SYMBOL"]],
                             pct_cor = pct_cor,
                             tpm_cor = tpm_cor)

df_long <- gather(correlation_df, key=Correlation, value= Spearman_Cor, c(2,3)) %>% 
  mutate(Gene = factor(Gene, levels = unique(Gene)))

ggplot(df_long, aes(x = Gene, y = Spearman_Cor, color = Correlation, group = Correlation)) + geom_line() + geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=c(0.7, -0.7), linetype="dashed", color = "green")

## Could also add counts here for number that reach the thresholds
correlation_df[[3]] > 0.7

iso_pct <- t(iso_pct) %>% as.data.frame() %>% clean_names()
iso_tpm <- t(iso_tpm) %>% as.data.frame() %>% clean_names()

### Show which %Iso are clustered, then remove those which are highly correlated (UMAP/Hierachacal clustering)
corclust_pct <- corclust(iso_pct[,phen_exp %>% dplyr::select(ends_with("_ci")) %>% colnames()] %>% na.omit(), cl = NULL, method = "complete")
plot(corclust_pct, mincor=0.7) 
clusters0.7_pct <- cvtree(corclust_pct, mincor = 0.7)
iso_removed <- as_tibble(clusters0.7_pct$correlations, rownames="Isoforms")[as_tibble(clusters0.7_pct$correlations, rownames="Isoforms")[2] %>% duplicated(),1]

#corclust_tpm <- corclust(iso_tpm[,phen_exp %>% dplyr::select(ends_with("_ci")) %>% colnames()] %>% na.omit(), cl = NULL, method = "complete")
#plot(corclust_tpm, mincor=0.7) 
#This is code for TPM if needs be

phen_exp <- phen_exp[,!colnames(phen_exp) %in% c(iso_removed[[1]], sub("\\_.*", "", iso_removed[[1]]))]

### Add in summary of number of variables removed from each pre-processing step

##### (3) Models

#################
### Model 1 Predicted early/late stage breast cancer  
#################

data1 <- phen_exp[-which(is.na(phen_exp$stage_binary)),] 

set.seed(3832)
trainIndex1 <- data1$stage_binary %>% createDataPartition(list=F, p=.8)

train <- data1[ trainIndex1,] %>% dplyr::select(-stage, -os, -os_time, -pam50) %>% na.omit()
test  <- data1[-trainIndex1,] %>% dplyr::select(-stage, -os, -os_time, -pam50) %>% na.omit()


logistic1 <- glm(stage_binary ~., 
                 data = train, 
                 family = "binomial")

logistic2 <- glm(stage_binary ~., 
                 data = train %>% dplyr::select(ends_with("_ci"), stage_binary),
                 family = "binomial")

sigiso <- filter(summary(logistic2)$coefficients %>% as.data.frame(),  `Pr(>|z|)` < 0.05)  %>% rownames()

logistic3 <- glm(stage_binary ~., 
                 data = train %>% dplyr::select(all_of(sigiso), stage_binary),
                 family = "binomial")

summary(logistic3)

roc(train$stage_binary ~ predict(logistic1, train, type = "response"),
    plot=T, 
    print.auc = T,
    main = "Model 1 Train")

roc(test$stage_binary ~ predict(logistic1, test, type = "response"),
    plot=T, 
    print.auc = T,
    main = "Model 1 Test")

roc(train$stage_binary ~ predict(logistic2, train, type = "response"),
    plot=T, 
    print.auc = T,
    main = "Model 2 Train")

roc(test$stage_binary ~ predict(logistic2, test, type = "response"),
    plot=T, 
    print.auc = T,
    main = "Model 2 Test")

roc <- roc(train$stage_binary ~ predict(logistic3, train, type = "response"),
    plot=T, 
    print.auc = T,
    main = "Model 3 Train")

roc(test$stage_binary ~ predict(logistic3, test, type = "response"),
    plot=T, 
    print.auc = T,
    main = "Model 3 Test") 

threshold <- coords(roc, "best", ret = "threshold")
confusionMatrix(as.factor(test$stage_binary), as.factor(ifelse(predict(logistic3, test, type = "response") < threshold[[1]], 0, 1)))





###xgboost models

tune_control <- trainControl("cv", 
                             number = 5,
                             verboseIter = FALSE,
                             allowParallel = TRUE)

final_grid <- expand.grid(
  nrounds = 100,
  eta = 0.01,
  max_depth = 2,
  gamma = 1,
  colsample_bytree = 0.5,
  min_child_weight = 0,
  subsample = 1
)


xgbmodel1 <- train(as.factor(stage_binary) ~.,
     data = train %>% dplyr::select(ends_with("_ci"), stage_binary),
     method = "xgbTree",
     trControl = tune_control,
     tuneGrid = final_grid,
     metric = "Kappa")


roc(train$stage_binary ~ predict(xgbmodel1$finalModel, train %>% dplyr::select(xgbmodel1[["finalModel"]][["feature_names"]]) %>% as.matrix(), type="response"),
    plot=T, 
    print.auc = T)
roc(test$stage_binary ~ predict(xgbmodel1$finalModel, test %>% dplyr::select(xgbmodel1[["finalModel"]][["feature_names"]]) %>% as.matrix(), type="response"),
    plot=T, 
    print.auc = T)

shap_values <- shap.values(xgb_model = xgbmodel1$finalModel, X_train = train %>% dplyr::select(xgbmodel1[["finalModel"]][["feature_names"]]) %>% as.matrix())
shap_long <- shap.prep(xgb_model = xgbmodel1$finalModel, X_train = train %>% dplyr::select(xgbmodel1[["finalModel"]][["feature_names"]]) %>% as.matrix())
shap.plot.summary.wrap1(xgbmodel1$finalModel, train %>% dplyr::select(xgbmodel1[["finalModel"]][["feature_names"]]) %>% as.matrix(), top_n=20)

shap.plot.dependence(data_long = shap_long, "icoslg_ci")
shap.plot.dependence(data_long = shap_long, "cct7_ci")
shap.plot.dependence(data_long = shap_long, "prr3_ci")

selectediso_data <- data.frame(stage_binary = data1$stage_binary, icoslg_ci = data1$icoslg_ci,
                               cct7_ci = data1$cct7_ci, prr3_ci = data1$prr3_ci)
selectediso_data$icoslg_ci_bin <- ifelse(selectediso_data$icoslg_ci < 20, 0, 1)
selectediso_data$prr3_ci_bin <- ifelse(selectediso_data$prr3_ci < 15, 0, 1)

train <- selectediso_data[ trainIndex1,] %>% na.omit()
test  <- selectediso_data[-trainIndex1,] %>% na.omit()

 
logistic4 <- glm(stage_binary ~ icoslg_ci, 
                   data = train,
                   family = "binomial")

logistic5 <- glm(stage_binary ~ prr3_ci_bin, 
                 data = train,
                 family = "binomial")

logistic6 <- glm(stage_binary ~ cct7_ci, 
                 data = train,
                 family = "binomial")


for (i in 4)
  
roc <- roc(train$stage_binary ~ predict(logistic4, train, type = "response"),
    plot=T, 
    print.auc = T)
threshold <- coords(roc, "best", ret = "threshold")
confusionMatrix(as.factor(test$stage_binary), as.factor(ifelse(predict(logistic4, test, type = "response") < threshold[[1]], 0, 1))) %>% print()

roc <- roc(train$stage_binary ~ predict(logistic5, train, type = "response"),
           plot=T, 
           print.auc = T)
threshold <- coords(roc, "best", ret = "threshold")
confusionMatrix(as.factor(test$stage_binary), as.factor(ifelse(predict(logistic5, test, type = "response") < threshold[[1]], 0, 1))) %>% print()

roc <- roc(train$stage_binary ~ predict(logistic6, train, type = "response"),
           plot=T, 
           print.auc = T)
threshold <- coords(roc, "best", ret = "threshold")
confusionMatrix(as.factor(test$stage_binary), as.factor(ifelse(predict(logistic6, test, type = "response") < threshold[[1]], 0, 1))) %>% print()


### Now for subtype
data2 <- phen_exp[-which(is.na(phen_exp$pam50)),] 
data2 <- data2[data2$pam50 != "Normal",] 
data2$pam50 <- as.factor(data2$pam50)

set.seed(3832)
trainIndex2 <- data2$pam50 %>% createDataPartition(list=F, p=.8)

train <- data2[ trainIndex1,] %>% dplyr::select(-stage_binary, -os, -os_time) %>% na.omit()
test  <- data2[-trainIndex1,] %>% dplyr::select(-stage_binary, -os, -os_time) %>% na.omit()

multilog1 <- multinom(pam50~.,
                      data=train)
confusionMatrix(test$pam50, predict(multilog1, test))

multilog2 <- multinom(pam50~.,
                      data = train %>% dplyr::select(ends_with("_ci"), pam50, stage))
confusionMatrix(test$pam50, predict(multilog2, test))


z <- summary(multilog2)$coefficients/summary(multilog2)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
sigiso_subtype <- p[,colSums(p < 0.05) > 0] %>% as.data.frame() %>% dplyr::select(ends_with("_ci")) %>% colnames()

multilog3 <- multinom(pam50~.,
                      data = train %>% dplyr::select(sigiso_subtype, pam50))
confusionMatrix(test$pam50, predict(multilog3, test))

set.seed(3445)
xgbmodel <- train(pam50 ~.,
                  data = dplyr::select(train, ends_with("ci"), pam50),
                  method = "xgbTree",
                  trControl = trainControl("cv", number = 5))

confusionMatrix(predict(xgbmodel, test) , test$pam50)



##Shap_value manually
shap_contrib <- predict(xgbmodel$finalModel, dplyr::select(train, ends_with("ci")) %>% as.matrix(), predcontrib = TRUE)
shap_contrib <- as.data.frame(shap_contrib[[1]])
BIAS0 <- shap_contrib[, ncol(shap_contrib)][1]
shap_contrib[, `:=`(BIAS, NULL)]
imp <- colMeans(abs(shap_contrib))
mean_shap_score <- imp[order(imp, decreasing = T)]
shap_values <- list(shap_score = shap_contrib, mean_shap_score = mean_shap_score, 
     BIAS0 = BIAS0)


shap_values <- shap.values(xgb_model = xgbmodel$finalModel, X_train = dplyr::select(train, ends_with("ci")) %>% as.matrix())
shap_long <- shap.prep(xgb_model = xgbmodel$finalModel, X_train = train %>% dplyr::select(xgbmodel[["finalModel"]][["feature_names"]]) %>% as.matrix())
shap.plot.summary.wrap1(xgbmodel$finalModel, train %>% dplyr::select(xgbmodel[["finalModel"]][["feature_names"]]) %>% as.matrix(), top_n=20)




#### Finally prediction of survival
data3 <- phen_exp[-which(is.na(phen_exp$os)),] 

set.seed(3832)
trainIndex3 <- data3$os %>% createDataPartition(list=F, p=.8)
dummy <- dummyVars(" ~ .", data=data3)
data3 <- data.frame(predict(dummy, newdata=data3))
data3$os <- as.numeric(data3$os)

train <- data3[ trainIndex3,] %>% dplyr::select(-stage1, - stage2, -stage3, -stage4) %>% na.omit()
test  <- data3[-trainIndex3,] %>% dplyr::select(-stage1, - stage2, -stage3, -stage4) %>% na.omit()


cox1 <- coxph(Surv(os_time, os) ~ ., data = train)
cox2 <- coxph(Surv(os_time, os) ~ ., data = train %>% dplyr::select(ends_with("_ci"), stage_binary, os, os_time, starts_with("pam50")))

sigiso <- filter(summary(cox2)$coefficients %>% as.data.frame(),  `Pr(>|z|)` < 0.05)  %>% rownames()

cox3 <- coxph(Surv(os_time, os) ~ ., data = train %>% dplyr::select(sigiso, os, os_time, starts_with("pam50")))

concordance(cox1)
concordance(cox2)
concordance(cox3)

# Read in tuning parameter from AZ work
load("/Users/am17168/Library/CloudStorage/OneDrive-UniversityofBristol/AZ/Isoform_proj/tune_hyper.RData")


### xgboost model
label.train <- ifelse(train$os == 1, train$os_time, -train$os_time)
label.test <- ifelse(test$os == 1, test$os_time, -test$os_time)

Dmat.train <- train %>% dplyr::select(ends_with("_ci"), stage_binary, age_days, starts_with("pam50")) %>%
  as.matrix() %>% xgb.DMatrix(label = label.train)

Dmat.test <- test %>% dplyr::select(ends_with("_ci"), stage_binary, os, os_time, starts_with("pam50")) %>%
  as.matrix() %>% xgb.DMatrix(label = label.test)

xgb1 <- xgboost(data = Dmat.train,
                params = tune_param_iso[6:12] %>% as.list(), 
                nrounds = tune_param_iso[[1]], 
                maximize = F, 
                eval_metric = "cox-nloglik")

shap_values <- shap.values(xgb_model = xgb1, X_train = Dmat.train)
shap_long <- shap.prep(xgb_model = xgb1, X_train = train %>% dplyr::select(ends_with("_ci"), stage_binary, age_days, starts_with("pam50")) %>% as.matrix())
shap.plot.summary.wrap1(xgb1, train %>% dplyr::select(ends_with("_ci"), stage_binary, age_days, starts_with("pam50")) %>% as.matrix(), top_n=20)

shap.plot.dependence(data_long = shap_long, "cct7_ci")
shap.plot.dependence(data_long = shap_long, "ccdc152_ci")
shap.plot.dependence(data_long = shap_long, "pnpla7_ci")
shap.plot.dependence(data_long = shap_long, "hoxc6_ci")



#Selected data
selectediso_data <- data3[, c("os", "os_time", "cct7_ci", "ccdc152_ci", "pnpla7_ci", "hoxc6_ci")]

selectediso_data$cct7_ci_bin <- ifelse(selectediso_data$cct7_ci <= 1, 0, 1)
selectediso_data$ccdc152_ci_bin <- ifelse(selectediso_data$ccdc152_ci < 60, 0, 1)
selectediso_data$pnpla7_ci_bin <- ifelse(selectediso_data$pnpla7_ci == 0, 0, 1)

train <- selectediso_data[ trainIndex3,] %>% na.omit()
test  <- selectediso_data[-trainIndex3,] %>% na.omit()

ggsurvplot(survfit(Surv(os_time, os)~cct7_ci_bin, data=test), 
           conf.int=TRUE, pval=TRUE, risk.table=TRUE, legend.labs=c(" Less than or equal to 1%", " More than 1%"), 
           title = "Km plot of cancerous isoform CCT7 at levels less than and greater than 1% of overall parent gene expression")

ggsurvplot(survfit(Surv(os_time, os)~ccdc152_ci_bin, data=test), 
           conf.int=TRUE, pval=TRUE, risk.table=TRUE, legend.labs=c(" Less than or equal to 60%", " More than 60%"), 
           title = "Km plot of cancerous isoform CCDC152 at levels less than and greater than 60% of overall parent gene expression")

ggsurvplot(survfit(Surv(os_time, os)~pnpla7_ci_bin, data=test), 
           conf.int=TRUE, pval=TRUE, risk.table=TRUE, legend.labs=c(" Less than or equal to 0%", " More than 0%"), 
           title = "Km plot of cancerous isoform PNPLA7 at levels less than and greater than 0% of overall parent gene expression")

