#################
#Cox Regression
#################

library(survival)
library(xgboost)
library(dplyr)
library(caret)
library(SHAPforxgboost)

load("Results/clean_data_all.RData")

#Remove variables if only have isoform OR gene expression
missS <- filter(BRCA_match_cc, !(CancerIso %in% iso))[["SYMBOL"]]
missI <- filter(BRCA_match_cc, !(SYMBOL %in% gene))[["CancerIso"]]

#Linear model
res.cox <- coxph(Surv(OS.time, OS) ~ ., data = clean_data_all)
summary(res.cox) 

#xgboost
trainIndex <- createDataPartition(clean_data_all$OS, 
                                  list=F,
                                  p=.8)
train <- clean_data_all[ trainIndex,]
test  <- clean_data_all[-trainIndex,]

label.train <- ifelse(train$OS == 1, train$OS.time, -train$OS.time)
label.test <- ifelse(test$OS == 1, test$OS.time, -test$OS.time)

Dmat.test <- test %>% as.matrix() %>% xgb.DMatrix(label = label.test)


#First build a model with default parametes
params <- list(booster = "gbtree", 
               objective = "survival:cox", 
               eta=0.3, #This controls the learning rate, between 0-1
               gamma=0, #This is the regularisation parameter, start with 0 and look at CV error rate, if train error >>> test erro then increase gamma
               max_depth=3, #Maximum depth of tree, larger the depth the more complex model which means more likely to overfit
               min_child_weight=1, #Block feature interaction to prevent overfitting
               subsample=1, #Controls the number of observations supplied to a tree, 0-1
               colsample_bytree=1 #Controls the number of features supplied to a tree, 0-1
)

xgbcv <- xgb.cv(data = Dmat.train,
                params = params, #The default parameters we have listed previously
                nrounds = 500, #Maximum number of rounds
                nfold = 5, #Number of folds,
                showsd = T, #Shows sd of cross validation
                stratified =T, #Statified by values of outcome labels
                print_every_n = 10, #prints message each 10 round rather than every round
                early_stopping_rounds = 20, #If doesn't improve for 20 round then the code will stop running
                maximize = F, #We must specify this when specify an early stopping round
)

#The model stopped at nrounds=21, meaning the first iteration showed the lowest test_cox_nloglik_mean
#The output also suggested that test and train errors were comparible meaning we will leave gamma=0

xgb1 <- xgb.train(data = Dmat.train,
                  params = params, 
                  nrounds = 10, 
                  #watchlist = list(val=dtest,train=dtrain), this is used when we have split our data previously to test and train
                  #print.every.n = 10, default is 1 so when only doing 10 rounds this line of code is not necessary 
                  maximize = F, 
                  eval_metric = "cox-nloglik")
importance <- xgb.importance(feature_names = colnames(Dmat.train), model = xgb1)
xgb.plot.importance(importance_matrix = importance)
