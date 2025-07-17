### Simulation of multiple DAGs 
library(igraph)
DAGs = function(n,p,example,K,pi){
  DAGs = list()
  N = n * K
  if (example %in% c("hub", "hub_mini", "hub_unbalance1", "hub_unbalance2")) {
    for (k in 1:K) {
      true_Adjace = matrix(0, p, p)
      cut.point.p = pi
      
      if (example == "hub") {
        true_Adjace[1, 2:cut.point.p] = sample(c(-1, 1), (cut.point.p - 1), replace = TRUE) *
          runif(cut.point.p - 1, 0.5, 1.5)
      } else if (example == "hub_mini") {
        if (k < (K/2+1)) {
          true_Adjace[1, 2:cut.point.p] = sample(c(-1, 1), (cut.point.p - 1), replace = TRUE) *
            runif(cut.point.p - 1, 0.5, 1)
        } else {
          true_Adjace[1, 2:cut.point.p] = runif(cut.point.p - 1, -0.5, 0.5)
        }
      } else if (example %in% c("hub_unbalance1", "hub_unbalance2")) {
        if (k == 1) {
          true_Adjace[1, 2:cut.point.p] = sample(c(-1, 1), (cut.point.p - 1), replace = TRUE) *
            runif(cut.point.p - 1, 0.5, 1)
          n = ifelse(example == "hub_unbalance1", round(N / 4), round(N * 3 / 4))
        } else {
          true_Adjace[1, 2:cut.point.p] = runif(cut.point.p - 1, -0.5, 0.5)
          n = ifelse(example == "hub_unbalance1", round(N * 3 / 4), round(N / 4))
        }
      }
      noise=matrix(0,n,p)
      DAG_sample=matrix(0,n,p)
      noise[,1] = matrix(runif(n,-3,3),n,1)
      noise[,2:cut.point.p] = matrix(runif(n*(cut.point.p-1), -3,3), n, cut.point.p-1)
      noise[,(cut.point.p+1):p] = matrix(runif(n*(p-cut.point.p), -3,3), n, p-cut.point.p)
      DAG_sample[,1] = noise[,1]
      DAG_sample[,2:cut.point.p] = DAG_sample[,1] %*% t(true_Adjace[1,2:cut.point.p]) + noise[,2:cut.point.p]
      DAG_sample[,(cut.point.p+1):p] = noise[,(cut.point.p+1):p]
      DAGs[[k]] = list(DAG_sample=DAG_sample,true_Adjace=true_Adjace)
    }
  }
  else{
    startgraph = make_empty_graph(3)
    A = barabasi.game(n = p, m = pi, start.graph = startgraph)
    A = as.matrix(as_adjacency_matrix(A))
    A = t(A)
    w = which(A != 0)
    
    # Main loop
    for (k in 1:K) {
      if (example == "BA") {
        A[w] = runif(length(w), 0.5, 1.5) * sample(c(-1, 1), length(w), replace = TRUE) 
      } else if (example == "BA_mini") {
        if (k < (K/2+1)) {
          A[w] = sample(c(-1, 1), length(w), replace = TRUE) * runif(length(w), 0.5, 1)
        } else {
          A[w] = runif(length(w), -0.5, 0.5)
        } 
      } else if (example %in% c("BA_unbalance1", "BA_unbalance2")) {
        if (k == 1) {
          A[w] = sample(c(-1, 1), length(w), replace = TRUE) * runif(length(w), 0.5, 1)
          n = ifelse(example == "BA_unbalance1", round(N / 4), round(N * 3 / 4))
        } else {
          A[w] = runif(length(w), -0.5, 0.5)
          n = ifelse(example == "BA_unbalance1", round(N * 3 / 4), round(N / 4))
        } 
      }
      noise = matrix(runif(n * p, -3, 3), n, p)
      X = matrix(0, n, p)
      
      for (i in 1:p) {
        parent = which(A[, i] != 0)
        if (length(parent) == 0) {
          X[, i] = noise[, i]
        } else {
          X[, i] = apply(A[, i] * t(X), 2, sum) + noise[, i]
        }
      }
      
      DAGs[[k]] = list(DAG_sample = X, true_Adjace = A)
    }
  }
  
  return(DAGs)
}


gaussian.dag = function(n,p,K,pi,m){
  DAGs = list()
  for (k in 1:K) {
    true_Adjace = matrix(0, p, p)
    cut.point.p = pi
    true_Adjace[1, 2:cut.point.p] = sample(c(-1, 1), (cut.point.p - 1), replace = TRUE) *
      runif(cut.point.p - 1, 0.5, 1.5)
    noise=matrix(0,n,p)
    DAG_sample=matrix(0,n,p)
     
    
    num_uni = p - m
    uniform_cols <- sort(sample(1:p, num_uni))
    
    noise[, uniform_cols] <- matrix(runif(n*(num_uni), -3,3), n, num_uni)
    gaussian_cols <- setdiff(1:p, uniform_cols)
    noise[, gaussian_cols] <- matrix(rnorm(n * length(gaussian_cols), 0, 1),
                                   n, length(gaussian_cols))
    
    DAG_sample[,1] = noise[,1]
    DAG_sample[,2:cut.point.p] = DAG_sample[,1] %*% t(true_Adjace[1,2:cut.point.p]) + noise[,2:cut.point.p]
    DAG_sample[,(cut.point.p+1):p] = noise[,(cut.point.p+1):p]
    DAGs[[k]] = list(DAG_sample=DAG_sample,true_Adjace=true_Adjace)
  }
  return(DAGs)
}
