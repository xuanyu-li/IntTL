Compute.reg.matrix = function(dags,hardth,K,tuning){
  n=nrow(dags[[1]]$DAG_sample)
  p=ncol(dags[[1]]$DAG_sample)
  Reg.matrix = matrix(0,K*(p-1), p)
  idx_num = 10
  for (j in 1:p) {
    if (p>80){
      ### Preparing matrix and vectors and groups
      cors <- abs(cor(dags[[1]]$DAG_sample[,-j], dags[[1]]$DAG_sample[,j]))
      top_100_idx <- order(cors, decreasing = TRUE)[1:idx_num]
      group = rep(1:idx_num , each = K)
      x = dags[[1]]$DAG_sample[,top_100_idx]
      if (K == 1) {
        x = x
      }   else {
        for (i in 2:K) {
          x = rbind(cbind(x,matrix(0,nrow=nrow(x),ncol=idx_num)), cbind(matrix(0,nrow=n,ncol=ncol(x)),dags[[i]]$DAG_sample[,top_100_idx]))
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
        permute = seq(1,(K-1)*idx_num+1,idx_num)
        for (i in 2:idx_num) {
          permute = c(permute,seq(i,(K-1)*idx_num+i,idx_num))
        }
        
        x = x[,permute]
      }
      
      
      
      
      
      
      #### group lasso
      
      #  gglasso.result = cv.gglasso(x = x, y = y,group = group,loss = "ls",pred.loss = "L2",lambda = sqrt(log(max(p,n))/n)*c(0.1,0.4,0.7,1),nfolds = 3)
      # gglasso.result = gglasso(x = x, y = y,group = group,loss = "ls", lambda = sqrt(log(max(p,n))/n)*c(0.1,0.4,0.7,1))
      gglasso.result = gglasso(x = x, y = y,group = group,loss = "ls", lambda = tuning * sqrt(log(max(p,n))/n))
      
      ## threshold
      
      med = sapply(1:(idx_num), function(j) median(abs(gglasso.result$beta[1:K + K*(j-1)])))
      
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
      
      
      
      
      all_top_100_idx = top_100_idx
      if (K != 1) {
        for (k in 2:K){
          all_top_100_idx = c(all_top_100_idx, top_100_idx+p-1)
        }
      }
      ### permute back
    if ( K != 1) {
      Reg.matrix[,j][all_top_100_idx] = gglasso.beta[invPerm(permute)]
    } else {
      Reg.matrix[,j][all_top_100_idx] = gglasso.beta
    }
    }
   else{
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
 
Compute.dcov_p=function(index.Y,added.Reg.matrix,dags,K){
  ### Input:
  ### index.Y index of residual corresponding to response variable
  ### reg.coeff coefficient of regression, which is of length K*(p-1)
  ### comp.X n\times p sample matrix; p>=2
  replicates = 500
  p=ncol(dags[[1]]$DAG_sample)
  n=nrow(dags[[1]]$DAG_sample)
  if (p> 90) {
    cors <- abs(cor(dags[[1]]$DAG_sample[,-index.Y], dags[[1]]$DAG_sample[,index.Y]))
    top_100_idx <- order(cors, decreasing = TRUE)[1:90]
    # Subset the x matrix to those top 100 variables
    top_100_idx <- sort(top_100_idx)
    reg.coeff = added.Reg.matrix[,index.Y]
    comp.dcov.Kcol = matrix(0,90,K)
    for (i in 1:K) {
      comp.resid = dags[[i]]$DAG_sample %*% reg.coeff[1:p + p*(i-1)] 
      other.X = dags[[i]]$DAG_sample[,-index.Y]
      if (p>2){
        comp.dcov=sapply(top_100_idx, function(j) energy::dcov.test(comp.resid, other.X[,j], R=replicates)$p.value)
      } else{
        comp.dcov=energy::dcov.test(comp.resid,other.X, R=replicates)$p.value
      }
      comp.dcov.Kcol[,i] = comp.dcov
    }
  }
  else{
    reg.coeff = added.Reg.matrix[,index.Y]
    comp.dcov.Kcol = matrix(0,p-1,K)
    for (i in 1:K) {
      comp.resid = dags[[i]]$DAG_sample %*% reg.coeff[1:p + p*(i-1)] 
      other.X = dags[[i]]$DAG_sample[,-index.Y]
      if (p>2){
        comp.dcov=sapply(1:(p-1), function(j) energy::dcov.test(comp.resid, other.X[,j], R=replicates)$p.value)
      } else{
        comp.dcov=energy::dcov.test(comp.resid,other.X, R=replicates)$p.value
      }
      comp.dcov.Kcol[,i] = comp.dcov
    }
  }
  
  
  if (K != 1) {
    comp.dcov.Kcol = apply(comp.dcov.Kcol, 1, sum) / K
  }
  pvalue = comp.dcov.Kcol
  pvalue = mean(sort(pvalue)[1:ceiling(length(p)/100)])
  return(pvalue)
}
