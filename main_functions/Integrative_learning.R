### Integrative Learning
library(dcov)
library(Matrix)
library(gglasso)
library(permute)
 

Compute.reg.matrix = function(dags,hardth,K,tuning){
  n=nrow(dags[[1]]$DAG_sample)
  p=ncol(dags[[1]]$DAG_sample)
  Reg.matrix = matrix(0,K*(p-1), p)
  for (j in 1:p) {
    ### Preparing matrix and vectors and groups
    group = rep(1:(p - 1), each = K)
    x = dags[[1]]$DAG_sample[,-j]
    if (K == 1) {
      x = x
      if (is.vector(x)) {
        x = as.matrix(x)
      }
    }   else if (is.vector(x)) {
      x = as.matrix(x)
      for (i in 2:K) {
        x = rbind(cbind(x,matrix(0,nrow=nrow(x),ncol=ncol(as.matrix(dags[[i]]$DAG_sample[,-j])))),    cbind(matrix(0,nrow=nrow(as.matrix(dags[[i]]$DAG_sample[,-j])),ncol=ncol(x)),dags[[i]]$DAG_sample[,-j]))
      }
    } else {
      for (i in 2:K) {
        x = rbind(cbind(x,matrix(0,nrow=nrow(x),ncol=ncol(dags[[i]]$DAG_sample[,-j]))),    cbind(matrix(0,nrow=nrow(dags[[i]]$DAG_sample[,-j]),ncol=ncol(x)),dags[[i]]$DAG_sample[,-j]))
      }
    }
    
    y = dags[[1]]$DAG_sample[,j]
    if (K != 1) {
      for (i in 2:K) {
        y = c(y,dags[[i]]$DAG_sample[,j])
      }
    }
    
    ## Since this package only accept group variable non-decreasing, we permute the columns of X   
    if ( K != 1) {
      permute = seq(1,(K-1)*(p-1)+1,(p-1))
      if (p == 2) {
        permute = permute
      } else{
        for (i in 2:(p-1)) {
          permute = c(permute,seq(i,(K-1)*(p-1)+i,(p-1)) )
        }
      }
      x = x[,permute]
    }
    
    
    
    
    
    
    #### group lasso
    
    #  gglasso.result = cv.gglasso(x = x, y = y,group = group,loss = "ls",pred.loss = "L2",lambda = sqrt(log(max(p,n))/n)*c(0.1,0.4,0.7,1),nfolds = 3)
    # gglasso.result = gglasso(x = x, y = y,group = group,loss = "ls", lambda = sqrt(log(max(p,n))/n)*c(0.1,0.4,0.7,1))
    gglasso.result = gglasso(x = x, y = y,group = group,loss = "ls", lambda = tuning * sqrt(log(max(p,n))/n))
    
    ## threshold
    
    med = sapply(1:(p-1), function(j) median(abs(gglasso.result$beta[1:K + K*(j-1)])))
    
    index.g = K*which(med > hardth)
    index.gall = index.g
    if (K != 1) {
      for (k in 1:(K-1)) {
        index.gall = c(index.gall, index.g - k)
      }
    }
    
    gglasso.beta = rep(0,length(gglasso.result$beta))
    gglasso.beta[sort(index.gall)] = gglasso.result$beta[sort(index.gall)] 
    
    
    
    ## Adaptive lasso
    ome = med
    ome[med < hardth] = 0.000001
    gglasso.result = gglasso(x = x, y = y,group = group,loss = "ls", lambda =   tuning * sqrt(log(max(p,n))/n)/20, pf = ome^(-1))
    #gglasso.beta = gglasso.result$beta
    
    ## thresholding again
    med = sapply(1:(p-1), function(j) median(abs(gglasso.result$beta[1:K + K*(j-1)])))
    
    index.g = K*which(med > hardth)
    index.gall = index.g
    if (K != 1) {
      for (k in 1:(K-1)) {
        index.gall = c(index.gall, index.g - k)
      }
    }
    
    gglasso.beta = rep(0,length(gglasso.result$beta))
    gglasso.beta[sort(index.gall)] = gglasso.result$beta[sort(index.gall)] 
     
    
    
    
    ### permute back
    if ( K != 1) {
      Reg.matrix[,j] = gglasso.beta[invPerm(permute)]
    } else {
      Reg.matrix[,j] = gglasso.beta
    }
  }
  added.Reg.matrix = matrix(0,K*p, p)
  for (j in 1:p) {
    indexes = seq(j, K*p-p+j, p)
    added.Reg.matrix[indexes, j] = -1
    added.Reg.matrix[-indexes, j ] = Reg.matrix[,j]
  }
  
  
  
  #added.Reg.matrix[abs(added.Reg.matrix) < hardth] = 0
  return(-added.Reg.matrix)
}
  

rep.dcovtest = function(index.Y,added.Reg.matrix,dags,K,replicates){
  p=ncol(dags[[1]]$DAG_sample)
  n=nrow(dags[[1]]$DAG_sample)
  reg.coeff = added.Reg.matrix[,index.Y]
  comp.dcov.Kcol = matrix(0,(p-1),replicates)
  for (i in 1:K) {
    comp.resid = dags[[i]]$DAG_sample %*% reg.coeff[1:p + p*(i-1)]
    other.X = dags[[i]]$DAG_sample[,-index.Y]
    if (p>2){
      comp.dcov=sapply(1:(p-1),function(j) energy::dcov.test(comp.resid, other.X[,j],R = replicates)$replicates)
      comp.dcov = t(comp.dcov)
    } else{
      comp.dcov=energy::dcov.test(comp.resid, other.X,R = replicates)$replicates
    }
    comp.dcov.Kcol = comp.dcov.Kcol + comp.dcov
  }
  
  return(comp.dcov.Kcol)
}


Compute.dcov_p=function(index.Y,added.Reg.matrix,dags,K){
  ### Input:
  ### index.Y index of residual corresponding to response variable
  ### reg.coeff coefficient of regression, which is of length K*(p-1)
  ### comp.X n\times p sample matrix; p>=2
  replicates = 700
  p=ncol(dags[[1]]$DAG_sample)
  n=nrow(dags[[1]]$DAG_sample)
  reg.coeff = added.Reg.matrix[,index.Y]
  non.zeros = which(reg.coeff != 0)
  comp.dcov.Kcol = matrix(0,(p-1),K)
  for (i in 1:K) {
    comp.resid = dags[[i]]$DAG_sample %*% reg.coeff[1:p + p*(i-1)] 
    other.X = dags[[i]]$DAG_sample[,-index.Y]
    if (p>2){
      comp.dcov=sapply(1:(p-1),function(j) energy::dcov.test(comp.resid, other.X[,j])$statistic)
    } else{
      comp.dcov=energy::dcov.test(comp.resid,other.X)$statistic
    }
    comp.dcov.Kcol[,i] = comp.dcov
  }
  if (K != 1) {
    comp.dcov.Kcol = apply(comp.dcov.Kcol, 1, sum)
  }
  
  permu = rep.dcovtest(index.Y,added.Reg.matrix,dags,K,replicates)
  pvalue = apply(permu, 2, function(x) x > comp.dcov.Kcol)
  if (p == 2) {
    pvalue = matrix(pvalue,1,replicates)
  }
  pvalue = apply(pvalue, 1, sum)/replicates
  pvalue = mean(sort(pvalue)[1:ceiling(length(non.zeros)/20)])
  return(pvalue)
}


TLLiNGAM=function(K,dags,hardth,criti.val_c,tuning){
  ### Input:
  ### K number of datesets
  ### X n\times p sample matrix
  ### hardth hard threshold of regression
  ### criti.val critical value of independence test based on distance covariance
  ### criti.val_c the constant of c_0 as suggested below
  ### as suggested the critical value can be c_0 * 1/n^(1/2) * (1+ log max(p,n) / K)^{1/2}
  origindags = dags
  check=FALSE
  n=nrow(dags[[1]]$DAG_sample)
  p=ncol(dags[[1]]$DAG_sample)
  K=length(dags)
  A=matrix(0,p,2) # store the information of layer
  A[,1]=1:p
  colnames(A)=c("Node","Layer")
  B = list()
  for (k in 1:K) {
    B[[k]] = matrix(0,p,p)
  }
  
  S=1:p
  # Identify A_0 layer
  
  ### June 22nd
  ## My regression vector is a (p-1)*K dim vector 
  
  added.Reg.matrix = Compute.reg.matrix(dags,hardth,K,tuning)
  non_zero_cols <- which(colSums(added.Reg.matrix != 0) > K)
  if (length(non_zero_cols)==0){
    node_S_select=S
    A[node_S_select,2]=0
    check=TRUE
    return(list(A=A,B=B,check=check))
  }
  dcov_pvalue_temp = sapply(non_zero_cols,function(j) Compute.dcov_p(j,added.Reg.matrix,dags,K))
  dcov_pvalue <- replace(numeric(p) + 1, non_zero_cols, dcov_pvalue_temp)  
  
  node_S_select=S[dcov_pvalue > criti.val_c]
  # if ((1 %in% node_S_select)) {
  #   x = 1
  # }
  
  if (length(node_S_select) == p) {
    criti.val = 2 * criti.val_c
    node_S_select = S[dcov_pvalue > criti.val]
    if (length(node_S_select) == p) {
      criti.val = min(dcov_pvalue)
      node_S_select = S[dcov_pvalue > criti.val]
    }
  }
  
  len_S = length(S)
  if (length(node_S_select)==0){
    check=TRUE
    return(list(A=A,B=B,check=check))
  }
  S=S[-node_S_select]
  
  for (k in 1:K) {
    B[[k]][S,node_S_select] = -added.Reg.matrix[1:len_S + (k-1)*len_S,][S,node_S_select]
  }
  
  
  
  t=1
  
  while(length(S)>1){
    len_S=length(S)
    
    
    ##### 
    
    for (k in 1:K) {
      dags[[k]]$DAG_sample = origindags[[k]]$DAG_sample[,S]
    }
    
    added.Reg.matrix = Compute.reg.matrix(dags,hardth,K,tuning)
    non_zero_cols <- which(colSums(added.Reg.matrix != 0) > K)
    if (length(non_zero_cols)==0){
      node_S_select=S
      A[node_S_select,2]=t
      check=TRUE
      return(list(A=A,B=B,check=check))
    }
    dcov_pvalue_temp = sapply(non_zero_cols,function(j) Compute.dcov_p(j,added.Reg.matrix,dags,K))
    dcov_pvalue <- replace(numeric(len_S) + 1, non_zero_cols, dcov_pvalue_temp)   
    
    node_S_select=S[dcov_pvalue > criti.val_c]
    index_node_S_select=which(dcov_pvalue>criti.val_c)
    
    if (length(node_S_select) == len_S) {
      criti.val = 2 * criti.val_c
      node_S_select = S[dcov_pvalue > criti.val]
      index_node_S_select=which(dcov_pvalue>criti.val)
      if (length(node_S_select) == len_S) {
        criti.val = min(dcov_pvalue)
        node_S_select = S[dcov_pvalue > criti.val]
        index_node_S_select=which(dcov_pvalue>criti.val)
      }
    }
    
    if (length(node_S_select)==0){
      node_S_select=S
      A[node_S_select,2]=t
      check=TRUE
      return(list(A=A,B=B,check=check))
    }
    A[node_S_select,2]=t
    S=setdiff(S,node_S_select)
    index_S=setdiff(c(1:len_S),index_node_S_select)
    
    if ( length(node_S_select)!=len_S ) {
      for (k in 1:K) {
        B[[k]][S,node_S_select] = -added.Reg.matrix[1:len_S + (k-1)*len_S,][index_S,index_node_S_select]
      }
    }
    
    
    #   B[S,node_S_select]=stdPreciMat[index_S,index_node_S_select]
    t=t+1
  }
  if (length(S)==1){
    A[S,2]=t
  }
  sumB = matrix(0,p, p)
  for (k in 1:K) {
    sumB = sumB + abs(B[[k]])^2
  }
  index.B = (sqrt(sumB/K) < hardth)
  for (k in 1:K) {
    B[[k]][index.B]=0
  }
  
  #B=Parent.set.refit(X,B)
  return(list(A=A,B=B,check=check))
}

