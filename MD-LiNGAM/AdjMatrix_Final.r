########### calculate weighted adjacency matrix by parent set ##################
Get_AdjMatrix=function(bin_adj,sample_matrix){
  ### Input:
  ### bin_adj: binary adjacency matrix
  ### sample_matrix: sample matrix (n\time p)
  
  ### Output:
  ### check: False: get estimated weighted adjacency matrix
  ###        True: fail to get estimated weighted adjacency matrix: estimated weighted adjacency matrix is equal to binary adjacency matrix
  ### estimated_adjace: estimated (weighted) adjacency matrix
  p=ncol(sample_matrix)
  n=nrow(sample_matrix)
  check=FALSE
  
  estimated_adjace=matrix(0,p,p)
  
  parent_size=apply(bin_adj,2,sum)
  max_parent_size=max(parent_size)
  if(max_parent_size>=n){
    check=TRUE
    estimated_adjace=bin_adj
  } else{
    for (i in 1:p){
      parents=which(bin_adj[,i]!=0)
      if (length(parents)>=1){
        estimated_adjace[parents,i]=MASS::ginv(t(sample_matrix[,parents])%*%sample_matrix[,parents])%*%t(sample_matrix[,parents])%*%sample_matrix[,i]
      }
    }
  }
  return(list(estimated_adjace=estimated_adjace,check=check))
}

