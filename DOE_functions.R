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
    # formul<- as.formula(paste0(y,"~."))
    # cf1 <- Boruta::Boruta(formul, data = df, doTrace = 2)
    # find_vars <- cf1$finalDecision[which(cf1$finalDecision!="Rejected")]
    # Imp_vars <- names(find_vars)
    Imp_vars <- var.sel.r2vim(x, y, no.runs = 50)
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

#' Variable selection using recurrent relative variable importance (r2VIM).
#'
#' Generates several random forests using all variables and different random
#' number seeds. For each run, the importance score is divided by the (absolute)
#' minimal importance score (relative importance scores). Variables are selected
#' if the minimal relative importance score is >= factor.
#'
#' Note: This function is a reimplementation of the R package \code{RFVarSelGWAS}.
#'
#' @inheritParams wrapper.rf
#' @param no.runs number of random forests to be generated
#' @param factor minimal relative importance score for a variable to be selected
#'
#' @return List with the following components:
#'   \itemize{
#'   \item \code{info} data.frame
#'   with information for each variable
#'   \itemize{
#'   \item vim.run.x = original variable importance (VIM) in run x
#'   \item rel.vim.run.x = relative VIM in run x
#'   \item rel.vim.min = minimal relative VIM over all runs
#'   \item rel.vim.med = median relative VIM over all runs
#'   \item selected = variable has been selected
#'   }
#'   \item \code{var} vector of selected variables
#'   }
#'
#'  @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # select variables
#' res = var.sel.r2vim(x = data[, -1], y = data[, 1], no.runs = 5, factor = 1)
#' res$var
#'
#' @export

var.sel.r2vim <- function(x, y, no.runs = 10, factor = 1, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                          no.threads = 1, method = "ranger", type = "regression") {
  
  ## importance for each run
  imp.all = NULL
  for (r in 1:no.runs) {
    print(paste("run", r))
    rf = wrapper.rf(x = x, y = y,
                    ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop, no.threads = no.threads,
                    method = method, type = type)
    imp.all = cbind(imp.all, get.vim(rf))
  }
  
  ## factors
  min.global = min(imp.all)
  if (min.global >= 0) {
    stop("Global minimal importance score is not negative!")
  }
  no.neg.min = 0
  fac = matrix(nrow = nrow(imp.all), ncol = ncol(imp.all),
               dimnames = dimnames(imp.all))
  for (i in 1:ncol(imp.all)) {
    x = imp.all[,i]
    min = min(x)
    if (min >= 0) {
      no.neg.min = no.neg.min + 1
      fac[,i] = x / abs(min.global)
    } else {
      fac[, i] = x / abs(min)
    }
  }
  if (no.neg.min > 0) {
    print(paste(no.neg.min, "runs with no negative importance score!"))
  }
  fac.min = apply(fac, 1, min)
  fac.med = apply(fac, 1, median)
  
  ## select variables
  ind.sel = as.numeric(fac.min >= factor)
  
  ## info about variables
  info = data.frame(imp.all, fac, fac.min, fac.med, ind.sel)
  colnames(info) = c(paste("vim.run.", 1:no.runs, sep = ""),
                     paste("rel.vim.run.", 1:no.runs, sep = ""),
                     "rel.vim.min", "rel.vim.median", "selected")
  return(list(info = info, var = sort(rownames(info)[info$selected == 1])))
}

#' Wrapper function to call random forests function.
#'
#' Provides an interface to different parallel implementations of the random
#' forest algorithm. Currently, only the \code{ranger} package is
#' supported.
#'
#' @param x matrix or data.frame of predictor variables with variables in
#'   columns and samples in rows (Note: missing values are not allowed).
#' @param y vector with values of phenotype variable (Note: will be converted to factor if
#'   classification mode is used).
#' @param ntree number of trees.
#' @param mtry.prop proportion of variables that should be used at each split.
#' @param nodesize.prop proportion of minimal number of samples in terminal
#'   nodes.
#' @param no.threads number of threads used for parallel execution.
#' @param method implementation to be used ("ranger").
#' @param type mode of prediction ("regression", "classification" or "probability").
#' @param ... further arguments needed for \code{\link[relVarId]{holdout.rf}} function only.
#'
#' @return An object of class \code{\link[ranger]{ranger}}.
#'
#' @import methods stats
#'
#' @export
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # regression
#' wrapper.rf(x = data[, -1], y = data[, 1],
#'            type = "regression", method = "ranger")

wrapper.rf <- function(x, y, ntree = 100, mtry.prop = 0.2, nodesize.prop = 0.1, no.threads = 1,
                       method = "ranger", type = "regression", ...) {
  
  ## check data
  if (length(y) != nrow(x)) {
    stop("length of y and number of rows in x are different")
  }
  
  if (any(is.na(x))) {
    stop("missing values are not allowed")
  }
  
  if (type %in% c("probability", "regression") & (is.character(y) | is.factor(y))) {
    stop("only numeric y allowed for probability or regression mode")
  }
  
  ## set global parameters
  nodesize = floor(nodesize.prop * nrow(x))
  mtry = floor(mtry.prop * ncol(x))
  if (mtry == 0) mtry = 1
  
  if (type == "classification") {
    #    print("in classification")
    y = as.factor(y)
  }
  
  ## run RF
  if (method == "ranger") {
    if (type == "probability") {
      y = as.factor(y)
      prob = TRUE
    } else {
      prob = FALSE
    }
    
    rf = ranger::ranger(data = data.frame(y, x),
                        dependent.variable.name = "y",
                        probability = prob,
                        importance = "permutation", scale.permutation.importance = FALSE,
                        num.trees = ntree,
                        mtry = mtry,
                        min.node.size = nodesize,
                        num.threads = no.threads,
                        write.forest = TRUE,
                        ...)
  } else {
    stop(paste("method", method, "undefined. Use 'ranger'."))
  }
  
  return(rf)
}


#' Error calculation.
#'
#' Calculates errors by comparing predictions with the true values. For
#' regression and probability mode, it will give root mean squared error (rmse) and
#' pseudo R-squared (rsq). For classification mode, overall accuracy (acc), overall
#' error (err), Matthews correlation coefficient (mcc), sensitivity (sens) and
#' specificity (spec) are returned.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}
#' @param true vector with true value for each sample
#' @param test.set matrix or data.frame of predictor variables for test set with variables in
#'   columns and samples in rows (Note: missing values are not allowed)
#' @inheritParams wrapper.rf
#'
#' @return numeric vector with two elements for regression and probability estimation (rmse, rsq) and
#' five elements for classification (acc, err, mcc, sens, spec)
#'
#' @export
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # random forest
#' rf = wrapper.rf(x = data[, -1], y = data[, 1],
#'                 type = "regression")
#'
#' # error
#' calculate.error(rf = rf, true = data[, 1])

calculate.error <- function(rf, true, test.set = NULL) {
  
  if (is(rf, "ranger")) {
    if (!is.null(test.set)) {
      pred = predict(rf, data = test.set)$predictions
    } else {
      pred = rf$predictions
    }
    if (rf$treetype == "Probability estimation") {
      pred = pred[, 2]
    }
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  
  if ((is(rf, "randomForest") && rf$type == "classification") |
      (is(rf, "ranger") && rf$treetype == "Classification")) {
    conf.matrix = table(pred = pred, true = true)
    tp = conf.matrix[2, 2]
    tn = conf.matrix[1, 1]
    fn = conf.matrix[2, 1]
    fp = conf.matrix[1, 2]
    
    ## accuracy
    acc = (tp + tn) / sum(conf.matrix)
    
    ## Matthews correlation coefficient
    mcc = (tp * tn - fp * fn) /
      sqrt( (tp + fn) * (tn + fp) * (tp + fp) * (tn + fn))
    
    ## sensitivity
    sens = tp / (tp + fn)
    
    ## specificity
    spec = tn / (fp + tn)
    
    error = c(err = 1 - acc, acc = acc, mcc = mcc, sens = sens, spec = spec)
  } else {
    mse = sum((pred - true)^2, na.rm = TRUE) / sum(!is.na(pred))
    
    ## pseudo R-squared uses sum of squared differences divided by n instead of variance!
    v = sum((true - mean(true))^2) / length(true)
    rsq = 1 - mse/v
    error = c(rmse = sqrt(mse), rsq = rsq)
  }
  
  return(error)
}


#' Get variable importance.
#'
#' Extracts variable importance depending on class of random forest object.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}
#'
#' @return numeric vector with importance value for each variable (in original order)
#'
#' @export

get.vim <- function(rf) {
  if (is(rf, "ranger")) {
    vim = ranger::importance(rf)
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  return(vim)
}