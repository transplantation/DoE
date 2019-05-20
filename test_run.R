
source("DOE_functions.R")
## identifying test IDs from cleaned data
Index_test <- testdata_index(heart.df.cleaned$ID, 2019)

cat_vars <- setdiff(pool_char_clean, "year1")

## encoding (method = c("numeric", "factor"))
df_num <- lapply(dfs, function(x) encode_cat(x,cat_vars,"numeric"))
df_cat <- lapply(dfs, function(x) encode_cat(x,cat_vars,"factor"))
All_data <- c(df_num, df_cat)

## assign names for each component in the list
names(All_data) <- c(paste0("NUM_",1:8), paste0("CAT_",1:8))

## remove unused objects
rm(dfs,df_num,df_cat)
#save3




## variable selection FFS
impute_no <- 16  #from 1 to 16
features_FFS <- vector(mode="list", impute_no)




for (i in 1:impute_no){
  cl <- makeCluster(4, type="SOCK")
  features_FFS[[i]] <- parSapply(cl, 1:5, select_vars, All_data[[i]], "year1", "FFS", Index_test, 2090)
  stopCluster(cl)
}

## variable selection LASSO
for (i in 1:impute_no){
  cl <- makeCluster(4, type="SOCK")
  features_FFS[[i]] <- parSapply(cl, 1:5, select_vars, All_data[[i]], "year1", "LASSO", Index_test, 2090)
  stopCluster(cl)
}


#features_LASSO <- vector(mode="list", impute_no)
#features_LASSO[[1]] <- select_vars(1, All_data[[9]], "year1", "LASSO", Index_test, seed=2090)

## Modeling
I_matrix <- matrix(c(rep(1:5,16),rep(1:16,each=5)),ncol=2,byrow=F)
iteration <- 5 #nrow(I_matrix)
cl <- makeCluster(4, type="SOCK")
modeling_result <- parSapply(cl, 1:iteration, modeling, All_data, "year1", I_matrix, Index_test, features_FFS, "glm", sampling_method="up") 
stopCluster(cl)

## result <- modeling(41, All_data, "year1", I_matrix, Index_test, features_LASSO, "glm", sampling_method = "up")
Performance <- matrix(NA, ncol=iteration, nrow=4)
rownames(Performance) <- c("auc","sen","spec","accu")
colnames(Performance) <- paste("Case", 1:iteration, sep="")
P <- rep(list(NA), iteration)

for (i in 1:iteration){
  Performance[,i] <- as.vector(modeling_result[1,i]$Performance[,1])
  P[[i]] <- modeling_result[2,i]$Predicted
}
 

