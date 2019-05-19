## dummy_maker() performs One-hot encoding (creates dummy variables) for categorical variables in a given data object
## input_data: a data object
## char_var: a vector of all independent categorical variables in the data
## Return Values: a data object after One-hot encoding is applied
dummy_maker <- function(input_data,char_var){
  for (i in 1:ncol(input_data)){
    if(names(input_data[i]) %in% char_var){
      # Use createDummyvars() function to create dummy variables for each categorical variable
      # The definition of createDummyvars() can be found in this file
      temp <- createDummyvars(input_data[i])
      names(temp) <- paste(names(input_data[i]),levels(as.factor(input_data[,i])),sep="_")
      input_data <- cbind(input_data,temp)
      input_data[,ncol(input_data)]<-NULL}
  }
  # We removed the dependent dummy variable in each categorical variable
  input_data <- input_data[-which(names(input_data) %in% char_var)]
  return(input_data)
}


## createDummyvars() function creates dummy variables for a given categorical variable
## data0: a column (vector, categoric) in the data
## Return Values: a data object after one-hot encoding is applied to the given variable
createDummyvars <- function(data0){
  dd<-as.data.frame(table(data0))
  dum_data<-as.data.frame(matrix(0, ncol = nrow(dd), nrow = nrow(data0)))
  names(dum_data) <- dd[,1]
  for(i in 1:ncol(dum_data)){
    dum_data[which(data0==names(dum_data)[i]),i]<-1
    dum_data[i]<-as.factor(dum_data[,i])
  }
  return(dum_data)
}


## encode_cat() retuns the data after the assigned encoding method for
## categorical variables
## df: data; method = c("numeric", "factor")
encode_cat <- function(df, cat_vars, method){
  cat_index <- which(colnames(df)%in%cat_vars)
  if (method=="numeric"){
    df[,cat_index] <- apply(df[,cat_index], 2, function(x) as.numeric(as.factor(x)))
  }else if (method=="factor"){
    variables <- colnames(df)[cat_index]
    df <- dummy_maker(df, variables)
  }
  return(df)
}


## testdata_index() creates IDs for test data
## n: number of rows from the entire dataset
## seed: random seed used in the function, the default is 2019
testdata_index <- function(n, seed=2019){
  set.seed <- seed
  index <- list()
  index0 <- c()
  train_list <- c()
  test_list <- c()
  for (i in 1:5){
    index[[i]] <- sample(setdiff(1:n, index0), floor(n/5))
    index0 <- c(index0, index[[i]])
  }
  return(index)
} 

## select_vars() does the variable selection for the data
## three variable selection methods are available: FFS, LASSO, Random Forest
## ii: the index to find the iith training data (we have 5 disjoint holdout datasets)
## df: the whole dataset
## y: name of the response variable
## method: FFS or LASSO or RF
## ID: the ID for the test data set 
## seed: random seed used in the function, the default is 201905
select_vars <- function(ii, df, y, method, ID, seed=201905){
  set.seed <- seed
  df <- df[!df$ID%in%ID[[ii]],]
  df <- df[,colnames(df)!="ID"]
  X <- df[,colnames(df)!=y]
  if (method=="FFS"){
    X$y <- df[,y]
    disc <- "MDL"
    threshold <- 0.001
    attrs.nominal <- numeric()
    result_temp <- Biocomb::select.fast.filter(X, disc.method=disc, threshold=threshold, attrs.nominal=attrs.nominal)
    Imp_vars <- as.character(result_temp$Biomarker)
  }else if (method=="LASSO"){
    '%ni%' <- Negate('%in%')
    X <- data.matrix(apply(X, 2, as.numeric))
    glmnet1 <- glmnet::cv.glmnet(x=X,y=as.factor(df[,y]),type.measure='auc', family="binomial")
    co <- coef(glmnet1,s = "lambda.1se")
    co <- as.matrix(co[-1,])
    Imp_vars <- row.names(co)[which((co[,1]!=0))]
    Imp_vars <- Imp_vars[Imp_vars %ni% '(Intercept)']
  }else if (method=="RF"){
    library(party)
    library(varImp)
    formul<- as.formula(paste0(y,"~."))
    cf1 <- Boruta::Boruta(formul, data = df, doTrace = 2)
    find_vars <- cf1$finalDecision[which(cf1$finalDecision!="Rejected")]
    Imp_vars <- names(find_vars)
  }
  return(list(Imp_vars))
}


## modeling() utilitizes a machine learning model in the caret package using cv with 5 folds
## The available models are glm, nnet, svmRadial, rf, gbm, earth, rpart, xgbTree, naive_bayes,
## treebag, lda, ranger
## it also evaluates the performance measures: AUC, Sensitivity, Specificity, Accuracy on the holdout set
## This function is used with parSapply() function in the package: snow
## All functions used in this function should be imported locally since multi-cores are used
## ii: the ii-th data indicator for different row in Index_matrix
## All_df: a list object contains all data with difference scenarios (16 components)
## TARGET: name of the response / dependent variable in the data
## Index_matrix: a matrix I(i,j) that helps to indicate which scenairo (j) and the corresponding 
## train and holdout data (i) and features used (j,i)
## assigned_seed: set a random seed
## methods_input="glm" : default
## Return Values: a list contains the following:
## 1. Performance measures from the holdout object
## 2. Predicted Survival Probabilities for the patients in the holdout object


modeling <- function(ii,All_df,TARGET,Index_matrix,Index_test,features,methods_input="glm", assigned_seed=2019, sampling_method=NULL){
  library(caret)
  library(AUC)
  library(MASS)
  
  formul <- as.formula(paste0(TARGET,"~."))
  
  test_no <- Index_matrix[ii,1]
  impute_no <- Index_matrix[ii,2]
  df <- All_df[[impute_no]]
  index_t<- which(df$ID%in%Index_test[[test_no]])
  df <- df[,c(features[[impute_no]][[test_no]], TARGET)]
  hold_out <- df[index_t,]
  traindata <- df[-index_t,]
  traindata$ID <- NULL
  
  set.seed((assigned_seed+ii))

  # 1: survival; 0: death
  traindata[,TARGET] <- as.factor(ifelse(traindata[,TARGET]==0, "Death", "Survival"))
  hold_out[,TARGET] <- as.factor(ifelse(hold_out[,TARGET]==0, "Death", "Survival"))
  
  # I used 5 fold cross validation 
  control_setting <- trainControl(method = "cv", number=5, sampling=sampling_method , summaryFunction = twoClassSummary, search="random", classProbs = TRUE, selectionFunction="tolerance")
  
  if (methods_input%in% c("glm", "nnet", "svmRadial")){
    result_model <- train(formul, data=traindata, method=methods_input, family="binomial",
                          trControl = control_setting, metric="ROC")
  }else if (methods_input%in%c("rf", "gbm", "earth", "rpart", "xgbTree", "naive_bayes", "ranger")){
    result_model <- train(formul,  data=traindata, method=methods_input,
                          trControl = control_setting, tuneLength=10, metric="ROC")
  }else if(methods_input=="lda"){
    result_model <- train(formul, data=traindata, method=methods_input, preProcess="pca", preProcOptions = list(method="BoxCox"),
                          trControl = control_setting, metric="ROC")
  }else if(methods_input=="treebag"){
    result_model <- train(formul, data=traindata, method=methods_input, family="binomial",
                          trControl = control_setting, tuneLength=10, metric="ROC")
  }
  
  resul_raw <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(hold_out)))
  colnames(resul_raw) <- c("TARGET", methods_input, "Probability")
  resul_raw$TARGET <- hold_out[TARGET]
  
  train_raw <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(traindata)))
  colnames(train_raw) <- c("TARGET", methods_input, "Probability")
  train_raw$TARGET <- traindata[TARGET]
  
  resul_pred_perf<-as.data.frame(matrix(NA, ncol = 1, nrow = 4))
  colnames(resul_pred_perf)<-c(methods_input)
  rownames(resul_pred_perf)<-c("auc","sen","spec","accu")
  train_auc <- NA
  
  #if(class(try(varImp(result_model),silent = TRUE))!="try-error"){
  train_raw$Probability <- predict(result_model, newdata=traindata, type="prob")[,2]
  train_raw[methods_input] <- predict(result_model, newdata=traindata, type="raw")
  #train_auc <- AUC::auc(roc(train_raw$Probability, traindata[,TARGET]))
  resul_raw[methods_input] <- predict(result_model, newdata=hold_out, type="raw")
  resul_raw$Probability <- predict(result_model, newdata=hold_out, type="prob")[,2]
  resul_pred_perf[1,1] <- AUC::auc(roc(resul_raw$Probability,hold_out[,TARGET]))
  resul_pred_perf[2,1] <- caret::sensitivity(resul_raw[,methods_input],hold_out[,TARGET])
  resul_pred_perf[3,1] <- caret::specificity(resul_raw[,methods_input],hold_out[,TARGET])
  resul_pred_perf[4,1] <- (as.data.frame(confusionMatrix(resul_raw[,methods_input], hold_out[,TARGET])$overall))[1,]
  #}
  return(list(Performance=resul_pred_perf, Predicted=resul_raw))
}



 