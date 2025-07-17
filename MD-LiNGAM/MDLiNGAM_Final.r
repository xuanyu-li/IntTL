######################################### MD-LiNGAM ############################
##############################################################
### estimate the order and parent set by package: highDLingam
### Copyright in Y.Samuel Wang (2020 Biometrika)
##############################################################
# library(doParallel)
# library(foreach)

# library(highDLingam)

# findGraphMulti
# output: $topOrder estimated toplogical ordering of variable: (the first is root
#         node and the last is sink node)
#         $parents: parents[[the order of node]]: .nodes
#MDLiNGAM=function(X,maxInDegree=3,degree=3,cutOffScaling=0.5){
MDLiNGAM=function(X,maxInDegree=3,degree=3,cutOffScaling=1){
  output=highDLingam::findGraphMulti(Y=X,maxInDegree=maxInDegree,degree=degree,cutOffScaling=cutOffScaling,verbose=F)
  ordered=output$topOrder
  parentSet=output$parents
  # print(ordered)
  # print(parentSet)
  
  p=ncol(X)
  # parentSet[[p]]=setdiff(parentSet[[p]],ordered[p])
  parentSet[[p]]=intersect(parentSet[[p]],ordered[1:(p-1)])
  estimated_adjace=matrix(0,p,p)
  
  # print(ordered)
  # print(parentSet)
  
  # use simple regression to esimate the connection strength(adjacency matrix) of DAG
  for(i in 2:p){
    response=ordered[i]
    regressor=parentSet[[i]]
    if (length(regressor)!=0){
      #model_lm=summary(lm(X[,response]~X[,regressor]))
      
      #sign_regress=model_lm$coef[-1,4]>0.05
      #store_regress=model_lm$coef[-1,1]
      #store_regress[sign_regress]=0
      
      #estimated_adjace[regressor,response]=store_regress
      estimated_adjace[regressor,response]=1
    }
  }
  return(list(estimated_adjace=estimated_adjace))
}

# dat=Hub_DAG(500,20,66,2)
# Y=dat$Hub_sample
# ncores=4
# cl=makeCluster(ncores)
# registerDoParallel(cl)
# MDLiNGAM(X=Y,maxInDegree=3,degree=4,cutOffScaling=0.8)
# stopCluster(cl)

