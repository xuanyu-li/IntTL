## running experiments
library("EvaluationMeasures")
library("huge")
library(parallel) 
source("DAGs_generate.R")
source("Integrative_learning.R")
source("CV_selection.R")
source("evaluation.R")
library(expm)
library(doParallel)
library(xtable)
library(JGL)
set.seed(1)
## This is the code for the table 1, which is the setting of n = 100, p = 50, K = 2
## Also, this is the strong signal case.
## So, run.exaple is always "hub" and “BA"





running<-function(iter,run.example,run.n,run.p,run.K, run.hardth,run.tuning,run.pi,criti.val_c){
  #set.seed(run.seed)
  dags = DAGs(run.n, run.p,run.example,run.K, run.pi)
  true_adjacent = list()
  for (k in 1:run.K) {
    true_adjacent[[k]] = dags[[k]]$true_Adjace
  } 
  start_tl=proc.time()
  estimated_adjace_tl=TLLiNGAM(run.K,dags, run.hardth, criti.val_c,run.tuning)$B
  diff_tl=proc.time()-start_tl
  store_tl=Evaluation.DAG(estimated_adjace_tl,true_adjacent,run.K)$Eval_result
  store_tllingam=c(store_tl)
  return(store_tllingam)
}

running_cv<-function(iter,run.example,run.n,run.p,run.K, run.hardth,run.tuning_grid,run.pi,criti.val_c, num_fold){
  #set.seed(run.seed)
  dags = DAGs(run.n, run.p,run.example,run.K, run.pi)
  true_adjacent = list()
  for (k in 1:run.K) {
    true_adjacent[[k]] = dags[[k]]$true_Adjace
  } 
  start_tl=proc.time()
  estimated_adjace_tl = cv_sect_TLLinGAM(run.tuning_grid, run.K, dags,run.hardth,criti.val_c,num_fold)$B
  diff_tl=proc.time()-start_tl
  store_tl=Evaluation.DAG(estimated_adjace_tl,true_adjacent,run.K)$Eval_result
  store_tllingam=c(store_tl)
  return(store_tllingam)
}

naive.running<-function(iter,run.example,run.n,run.p,run.K, run.hardth,run.tuning,run.pi,criti.val_c){
  dags = DAGs(run.n, run.p,run.example,run.K, run.pi)
  true_adjacent = list()
  for (k in 1:run.K) {
    true_adjacent[[k]] = dags[[k]]$true_Adjace
  }
  for (k in 2:run.K) {
    dags[[1]]$DAG_sample =  rbind(dags[[1]]$DAG_sample, dags[[k]]$DAG_sample)
  }
  dags = list(dags[[1]])
  ## The mulplication of run.K is because when we set run.tuning we devided run.K for uniformity, however, we only use one dataset. And now we correct this run.K coefficient.
  estimated_adjace_tl=TLLiNGAM(1,dags, run.hardth, criti.val_c,run.K*run.tuning)$B
  store_tl=Evaluation.DAG(estimated_adjace_tl,true_adjacent,1)$Eval_result
  store_tllingam=c(store_tl)
  return(store_tllingam)
}

running.GGL <- function(iter,run.example,run.n,run.p,run.K, run.hardth,run.tuning,run.pi,criti.val_c){ 
  dags = DAGs(run.n, run.p,run.example,run.K, run.pi)
  true_adjacent = list()
  for (k in 1:run.K) {
    mat = dags[[k]]$true_Adjace
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    true_adjacent[[k]] = mat
  }
  DAG.sample = list()
  for (k in 1:run.K) {
    DAG.sample[[k]] = dags[[k]]$DAG_sample
  }
  
  start_tl=proc.time()
  ggl.results = JGL(Y=DAG.sample,penalty="group",lambda1=0,lambda2=2*run.tuning/run.K *sqrt(log(max(run.p,run.n))/run.n),return.whole.theta=TRUE)
  estimated_adjace_tl = ggl.results$theta
  diff_tl=proc.time()-start_tl
  store_tl=Evaluation.DAG(estimated_adjace_tl,true_adjacent,run.K)$Eval_result
  store_tllingam=c(store_tl) 
  return(store_tllingam)
}

run_simu = function(run.n, run.p, run.K, run.example, run.pi,run.hardth, run.tuning,criti.val_c, num_fold, run.tuning_grid, cv_bool){
  iter = 50 
  clnum<-detectCores()
  clnum = 50
  start_time <- proc.time() 
  cl <- makeCluster(getOption("cl.cores", clnum))
  clusterEvalQ(cl,library(Matrix))
  clusterEvalQ(cl,library(gglasso))
  clusterEvalQ(cl,library(igraph))
  
  clusterExport(cl = cl, varlist = c('iter','running','DAGs', 'TLLiNGAM','Compute.reg.matrix','Compute.dcov_p','Evaluation.DAG','hammingDistance','rep.dcovtest', 'create_kfolds', 'kfold', 'cv_evaluate', 'cv_select', 'cv_sect_TLLinGAM'))
  if (cv_bool == TRUE){
    result = parSapply(cl,1:iter, running_cv, run.example = run.example, run.n = run.n, run.p = run.p, run.K = run.K ,run.hardth = run.hardth, run.tuning_grid = run.tuning_grid, run.pi = run.pi,criti.val_c = criti.val_c, num_fold = num_fold)
  }
  else{
    result = parSapply(cl,1:iter, running, run.example = run.example , run.n = run.n, run.p = run.p, run.K = run.K ,run.hardth = run.hardth, run.tuning = run.tuning, run.pi = run.pi,criti.val_c = criti.val_c)
  }
  
  stopCluster(cl)
  end_time <- proc.time()
  elapsed_time <- as.numeric((proc.time() - start_time)['elapsed'])
  
  result.mean = apply(result, 1, mean)
  result.mean = round(result.mean,3)
  
  result.std  = apply(result, 1, function(x) sqrt(var(x)))
  result.std = round(result.std,2)
  result.meanstd = paste(result.mean,result.std,sep = "(")
  result.meanstd = paste(result.meanstd,c(),sep = ")")
  result.meanstd = append(result.meanstd, "IntTL",0)
  result.meanstd = c(result.meanstd, elapsed_time/iter)
  col.names = c("Method","TPR","FDR","HM","MCC","time")
  
  data.result = data.frame(as.matrix(t(result.meanstd)) )
  colnames(data.result) = col.names
  
  ##  TL
  start_time <- proc.time() 
  result = c()
  cl <- makeCluster(getOption("cl.cores", clnum))
  clusterEvalQ(cl,library(Matrix))
  clusterEvalQ(cl,library(gglasso))
  clusterEvalQ(cl,library(igraph))
  clusterExport(cl = cl, varlist = c('iter','running','DAGs', 'TLLiNGAM','Compute.reg.matrix','Compute.dcov_p','Evaluation.DAG','hammingDistance','rep.dcovtest', 'create_kfolds', 'kfold', 'cv_evaluate', 'cv_select', 'cv_sect_TLLinGAM'))
  
  result = cbind(result,parSapply(cl,1:iter, running, run.example = run.example, run.n = run.n, run.p = run.p, run.K = 1,run.hardth = run.hardth, run.tuning = run.K*run.tuning, run.pi = run.pi,criti.val_c = criti.val_c))
  stopCluster(cl)
  end_time <- proc.time()
  elapsed_time <- as.numeric((proc.time() - start_time)['elapsed'])
  
  result.mean = apply(result, 1, mean)
  result.mean = round(result.mean,3)
  
  result.std  = apply(result, 1, function(x) sqrt(var(x)))
  result.std = round(result.std,2)
  result.meanstd = paste(result.mean,result.std,sep = "(")
  result.meanstd = paste(result.meanstd,c(),sep = ")")
  result.meanstd = append(result.meanstd, "TL",0)
  result.meanstd = c(result.meanstd, elapsed_time/iter)
  
  data.result = rbind(data.result, result.meanstd) 
  
  
  ## naive-TL
  start_time <- proc.time() 
  result = c()
  cl <- makeCluster(getOption("cl.cores", clnum))
  clusterEvalQ(cl,library(Matrix))
  clusterEvalQ(cl,library(gglasso))
  clusterEvalQ(cl,library(igraph))
  clusterExport(cl = cl, varlist = c('iter','running','DAGs', 'TLLiNGAM','Compute.reg.matrix','Compute.dcov_p','Evaluation.DAG','hammingDistance','rep.dcovtest', 'create_kfolds', 'kfold', 'cv_evaluate', 'cv_select', 'cv_sect_TLLinGAM'))
  result = cbind(result,parSapply(cl,1:iter, naive.running, run.example = run.example, run.n = run.n, run.p = run.p, run.K = run.K,run.hardth = run.hardth, run.tuning = run.tuning, run.pi = run.pi,criti.val_c = criti.val_c))
  stopCluster(cl)
  end_time <- proc.time()
  elapsed_time <- as.numeric((proc.time() - start_time)['elapsed'])
  
  result.mean = apply(result, 1, mean)
  result.mean = round(result.mean,3)
  
  result.std  = apply(result, 1, function(x) sqrt(var(x)))
  result.std = round(result.std,2)
  result.meanstd = paste(result.mean,result.std,sep = "(")
  result.meanstd = paste(result.meanstd,c(),sep = ")")
  result.meanstd = append(result.meanstd, "naive-TL",0)
  result.meanstd = c(result.meanstd, elapsed_time/iter)
  
  data.result = rbind(data.result, result.meanstd)
  ## GGL-moralized graph
  start_time <- proc.time() 
  result = c()
  cl <- makeCluster(getOption("cl.cores", clnum))
  clusterEvalQ(cl,library(Matrix))
  clusterEvalQ(cl,library(gglasso))
  clusterEvalQ(cl,library(igraph))
  clusterEvalQ(cl,library(JGL))
  clusterExport(cl = cl, varlist = c('iter','running','DAGs', 'TLLiNGAM','Compute.reg.matrix','Compute.dcov_p','Evaluation.DAG','hammingDistance','rep.dcovtest'))
  
  result = cbind(result,parSapply(cl,1:iter, running.GGL, run.example = run.example, run.n = run.n, run.p = run.p, run.K = 1,run.hardth = run.hardth, run.tuning = run.K*run.tuning, run.pi = run.pi,criti.val_c = criti.val_c))
  stopCluster(cl)
  end_time <- proc.time()
  elapsed_time <- as.numeric((proc.time() - start_time)['elapsed'])
  
  result.mean = apply(result, 1, mean)
  result.mean = round(result.mean,3)
  
  result.std  = apply(result, 1, function(x) sqrt(var(x)))
  result.std = round(result.std,2)
  result.meanstd = paste(result.mean,result.std,sep = "(")
  result.meanstd = paste(result.meanstd,c(),sep = ")")
  result.meanstd = append(result.meanstd, "GGL",0)
  result.meanstd = c(result.meanstd, elapsed_time/iter)
  
  data.result = rbind(data.result, result.meanstd) 
  
  xtable(data.result)
}






###  hub p 81 pi 2
run.n = 150
run.p = 81
run.K = 2
run.example = "hub"
run.pi = 8
run.hardth = 0.5*(log(run.p)/run.n/run.K)^(1/2)
run.tuning = 2
run.tuning_grid = c(0.2, 2, 20)
criti.val_c = 0.01
num_fold = 5
cv_bool = FALSE
run_simu(run.n, run.p, run.K, run.example, run.pi,run.hardth, run.tuning,criti.val_c, num_fold, run.tuning_grid,cv_bool)
 
###  hub p 81 pi 3
run.n = 150
run.p = 81
run.K = 2
run.example = "hub"
run.pi = 16
run.hardth = 0.5*(log(run.p)/run.n/run.K)^(1/2)
run.tuning = 2
run.tuning_grid = c(0.2, 2, 20)
criti.val_c = 0.01
num_fold = 5
cv_bool = FALSE
run_simu(run.n, run.p, run.K, run.example, run.pi,run.hardth, run.tuning,criti.val_c, num_fold, run.tuning_grid,cv_bool)

###  hub p 81 pi 4
run.n = 150
run.p = 81
run.K = 2
run.example = "hub"
run.pi = 24
run.hardth = 0.5*(log(run.p)/run.n/run.K)^(1/2)
run.tuning = 2
run.tuning_grid = c(0.2, 2, 20)
criti.val_c = 0.01
num_fold = 5
cv_bool = FALSE
run_simu(run.n, run.p, run.K, run.example, run.pi,run.hardth, run.tuning,criti.val_c, num_fold, run.tuning_grid,cv_bool)

source("Integrative_high.R")
###  hub p 243 pi 2
run.n = 200
run.p = 243
run.K = 2
run.example = "hub"
run.pi =  24
run.hardth = 0.5*(log(run.p)/run.n/run.K)^(1/2)
run.tuning = 2
run.tuning_grid = c(0.2, 2, 20)
criti.val_c = 0.01
num_fold = 5
cv_bool = FALSE
run_simu(run.n, run.p, run.K, run.example, run.pi,run.hardth, run.tuning,criti.val_c, num_fold, run.tuning_grid,cv_bool)

###  hub p 243 pi 3
run.n = 200
run.p = 243
run.K = 2
run.example = "hub"
run.pi =  48
run.hardth = 0.5*(log(run.p)/run.n/run.K)^(1/2)
run.tuning = 2
run.tuning_grid = c(0.2, 2, 20)
criti.val_c = 0.01
num_fold = 5
cv_bool = FALSE
run_simu(run.n, run.p, run.K, run.example, run.pi,run.hardth, run.tuning,criti.val_c, num_fold, run.tuning_grid,cv_bool)

###  hub p 243 pi 4
run.n = 200
run.p = 243
run.K = 2
run.example = "hub"
run.pi =  72
run.hardth = 0.5*(log(run.p)/run.n/run.K)^(1/2)
run.tuning = 2
run.tuning_grid = c(0.2, 2, 20)
criti.val_c = 0.01
num_fold = 5
cv_bool = FALSE
run_simu(run.n, run.p, run.K, run.example, run.pi,run.hardth, run.tuning,criti.val_c, num_fold, run.tuning_grid,cv_bool)