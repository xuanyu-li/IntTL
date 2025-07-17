# cross validate the tuning paramter selection
create_kfolds <- function(data, k = 5, seed = 123) {
  set.seed(seed)  # for reproducibility
  n <- nrow(data)
  indices <- sample(n)  # shuffle row indices
  folds <- split(indices, cut(seq_along(indices), breaks = k, labels = FALSE))
  return(folds)
}

kfold = function(dags,K, num_fold){
  dags_folds_train = list()
  dags_folds_test = list()
  
  folds = create_kfolds(dags[[1]]$DAG_sample, num_fold)
  for (j in 1:num_fold) {
    temp_fold_train = list()
    temp_fold_test = list()
    for (k in 1:K) {
      temp_fold_train[[k]] = list(DAG_sample = dags[[k]]$DAG_sample[-folds[[j]],])
      temp_fold_test[[k]] =  list(DAG_sample = dags[[k]]$DAG_sample[folds[[j]] ,])
    }
    dags_folds_train[[j]] = temp_fold_train
    dags_folds_test[[j]] = temp_fold_test
  }
  return(list(trainfolds = dags_folds_train, testfolds = dags_folds_test))
}

 

cv_evaluate = function(result, testfolds, num_fold, K){
  loss = numeric(num_fold)
  for (i in 1:num_fold) {
    estimated_B = result[,i]$B
    data=testfolds[[i]]
    for (k in 1:K) {
      loss[i] =  loss[i] + norm(data[[k]]$DAG_sample - data[[k]]$DAG_sample %*% estimated_B[[k]], type = "F")
    }
  }
  return(mean(loss))
}
  
  
cv_select = function(run.tuning_grid, K,dags,hardth,criti.val_c,num_fold){
  num_fold = 5
  evaluate = numeric(length(run.tuning_grid))
  for (run.tuning  in run.tuning_grid) {
    trainfolds = kfold(dags,K, num_fold)$trainfolds
    testfolds = kfold(dags,K, num_fold)$testfolds
    result = sapply(1:num_fold,function(j) TLLiNGAM(K,trainfolds[[j]],hardth,criti.val_c,run.tuning))
    evaluate[which(run.tuning_grid == run.tuning)] = cv_evaluate(result, testfolds, num_fold, K)
  }
  return(evaluate)
}

cv_sect_TLLinGAM = function(run.tuning_grid, K, dags,hardth,criti.val_c,num_fold){
  evaluate = cv_select(run.tuning_grid, K,dags,hardth,criti.val_c,num_fold)
  select_run = run.tuning_grid[which.min(evaluate)]
  result = TLLiNGAM(K, dags, hardth, criti.val_c,select_run)
  return(result)
}
