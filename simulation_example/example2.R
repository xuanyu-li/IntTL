# R script to run a simulation and produce a cross-fitting tuning parameter selection for the intTL method.

set.seed(2025)

library(Matrix)
library(igraph)
library(JGL)           # for the GGL baseline

# Correct sourcing paths and ensure functions are available
source("../simulation_example/DAGs_generate.R")
source("../main_functions/Integrative_learning.R")
source("../main_functions/CV_selection.R")
source("../main_functions/evaluation.R")
source("../MD-LiNGAM/MDLiNGAM_Final.r")
source("../MD-LiNGAM/AdjMatrix_Final.r")

#----------------------------
# Parameters for one illustrative run
#----------------------------
run.n <- 100      # sample size per dataset
run.p <- 27       # number of variables
run.K <- 2        # number of related datasets (integrative setting)
run.example <- "hub"
run.pi <- ceiling(2 * log(run.p))
run.hardth <- 0.5 * (log(run.p) / run.n / run.K)^(1/2)
run.tuning_grid = c(0.2, 2, 20)
criti.val_c <- 0.01

 

#----------------------------
# Cross-validation tuning selection
#----------------------------
run_cv <- function(run.n, run.p, run.K, run.example, run.pi, run.hardth, run.tuning_grid, criti.val_c, num_fold) {
  ## 1) generate K DAG datasets
  dags <- DAGs(run.n, run.p, run.example, run.K, run.pi)

  ## prepare true adjacency list for evaluation
  true_adjacent <- vector("list", run.K)
  for (k in seq_len(run.K)) true_adjacent[[k]] <- dags[[k]]$true_Adjace

  ## 2) Cross-validation for IntTL
  result <- cv_sect_TLLinGAM(run.tuning_grid, run.K, dags, run.hardth, criti.val_c, num_fold)
  eval_result <- Evaluation.DAG(result$B, true_adjacent, run.K)$Eval_result

  ## Return evaluation results
  return(eval_result)
}

#----------------------------
# Run cross-validation and show results
#----------------------------
cat("Running cross-validation for IntTL...\n")
cv_results <- run_cv(run.n, run.p, run.K, run.example, run.pi, run.hardth, run.tuning_grid, criti.val_c, num_fold = 5)
print(round(cv_results, 3))
