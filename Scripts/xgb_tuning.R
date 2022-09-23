param <- list(booster = "gbtree",
              objective = "survival:cox",
              max_depth = c(1,3,5),
              eta = c(.1, .3, 0.5),
              subsample = c(1, .7, .5),
              colsample_bytree = c(1, .7, 0.5),
              min_child_weight = c(0,1,2)
) %>% cross_df() %>% as.data.frame()

lowest_error_list <- list()

# Use randomly created parameters to create 10,000 XGBoost-models
for (row in 1:nrow(param)){
  xgbcv <- xgb.cv(data = Dmat.train,
                  params = param[row,] %>% as.list(), #The default parameters we have listed previously
                  nrounds = 500, #Maximum number of rounds
                  nfold = 5, #Number of folds,
                  showsd = T, #Shows sd of cross validation
                  stratified =T, #Statified by values of outcome labels
                  print_every_n = 10, #prints message each 10 round rather than every round
                  early_stopping_rounds = 20, #If doesn't improve for 20 round then the code will stop running
                  maximize = F, #We must specify this when specify an early stopping round
    )
  lowest_error <- as.data.frame(xgbcv$evaluation_log[xgbcv$best_iteration])
  lowest_error_list[[row]] <- lowest_error 
}

# Create object that contains all accuracy's
lowest_error_df = do.call(rbind, lowest_error_list)

# Bind columns of accuracy values and random hyperparameter values
randomsearch = cbind(lowest_error_df, param)

# Quickly display highest accuracy
tune_param_test <- randomsearch[which.min(randomsearch$train_cox_nloglik_mean),]

xgb1 <- xgboost(data = Dmat.train,
                params = tune_param[6:12] %>% as.list(), 
                nrounds = tune_param[[1]], 
                #watchlist = list(val=dtest,train=dtrain), this is used when we have split our data previously to test and train
                #print.every.n = 10, default is 1 so when only doing 10 rounds this line of code is not necessary 
                maximize = F, 
                eval_metric = "cox-nloglik")


### Now repeat but include gene expression as well
label.train <- ifelse(train$OS == 1, train$OS.time, -train$OS.time)
label.test <- ifelse(test$OS == 1, test$OS.time, -test$OS.time)
Dmat.train_gene <- train %>% as.matrix() %>% xgb.DMatrix(label = label.train)
Dmat.test_gene <- test %>% as.matrix() %>% xgb.DMatrix(label = label.test)

param <- list(booster = "gbtree",
              objective = "survival:cox",
              max_depth = c(1,3,5),
              eta = c(.1, .3, 0.5),
              subsample = c(1, .7, .5),
              colsample_bytree = c(1, .7, 0.5),
              min_child_weight = c(0,1,2)
) %>% cross_df() %>% as.data.frame()

lowest_error_list <- list()


n.cores <- 20
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
library(doParallel)

foreach(row = 1:nrow(param),
        .combine = 'rbind',
        .packages = c("xgboost", "dplyr")) %dopar% {
  xgbcv <- xgboost::xgb.cv(data = Dmat.train_gene,
                  params = as.list(param[row,]), #The default parameters we have listed previously
                  nrounds = 1500, #Maximum number of rounds
                  nfold = 5, #Number of folds,
                  showsd = T, #Shows sd of cross validation
                  stratified =T, #Statified by values of outcome labels
                  print_every_n = 10, #prints message each 10 round rather than every round
                  early_stopping_rounds = 20, #If doesn't improve for 20 round then the code will stop running
                  maximize = F, #We must specify this when specify an early stopping round
  ) %>% as.matrix()
  lowest_error <- as.data.frame(xgbcv$evaluation_log[xgbcv$best_iteration])
  lowest_error_list[[row]] <- lowest_error 
        }

parallel::stopCluster(cl = my.cluster)

lowest_error_list <- list()

for(row in 1:nrow(param)) {
          xgbcv <- xgboost::xgb.cv(data = Dmat.train_gene,
                                   params = as.list(param[row,]), #The default parameters we have listed previously
                                   nrounds = 1000, #Maximum number of rounds
                                   nfold = 5, #Number of folds,
                                   showsd = T, #Shows sd of cross validation
                                   stratified =T, #Statified by values of outcome labels
                                   print_every_n = 10, #prints message each 10 round rather than every round
                                   early_stopping_rounds = 20, #If doesn't improve for 20 round then the code will stop running
                                   maximize = F, #We must specify this when specify an early stopping round
          )
          lowest_error <- as.data.frame(xgbcv$evaluation_log[xgbcv$best_iteration])
          lowest_error_list[[row]] <- lowest_error 
}

# Create object that contains all accuracy's
lowest_error_df = do.call(rbind, lowest_error_list)

# Bind columns of accuracy values and random hyperparameter values
randomsearch = cbind(lowest_error_df, param)

# Quickly display highest accuracy
tune_param_gene <- randomsearch[which.min(randomsearch$train_cox_nloglik_mean),]
