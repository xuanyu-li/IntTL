## This script runs one data-generation and compares methods once.
## It is intentionally short and documented so new users can quickly run and inspect results.

## Usage: open this file in R and run source('example.R') 

set.seed(1)

library(Matrix)
library(igraph)
library(JGL)           # for the GGL baseline

source("DAGs_generate.R")                            # DAGs() (local to this folder)
source("../main_functions/Integrative_learning.R")  # TLLiNGAM() and related
source("../main_functions/CV_selection.R")          # CV utilities (kept for compatibility)
source("../main_functions/evaluation.R")            # Evaluation.DAG()


#----------------------------
# Parameters for one illustrative run
#----------------------------
run.n <- 100      # sample size per dataset
run.p <- 27       # number of variables
run.K <- 2        # number of related datasets (integrative setting)
run.example <- "hub"
run.pi <- ceiling(2 * log(run.p))
run.hardth <- 0.5 * (log(run.p) / run.n / run.K)^(1/2)
run.tuning <- 2
criti.val_c <- 0.01


#----------------------------
# Helper: run each method once and return results
#----------------------------
run_one <- function(run.n, run.p, run.K, run.example, run.pi, run.hardth, run.tuning, criti.val_c) {
  ## 1) generate K DAG datasets
  dags <- DAGs(run.n, run.p, run.example, run.K, run.pi)
  
  ## prepare true adjacency list for evaluation
  true_adjacent <- vector("list", run.K)
  for (k in seq_len(run.K)) true_adjacent[[k]] <- dags[[k]]$true_Adjace
  
  ## 2) IntTL (integrative TLLiNGAM) — run on all K datasets together
  est_inttl <- TLLiNGAM(run.K, dags, run.hardth, criti.val_c, run.tuning)$B
  eval_inttl <- Evaluation.DAG(est_inttl, true_adjacent, run.K)$Eval_result
  
  ## 3) TL (pooled single-dataset TLLiNGAM) — combine samples, run TLLiNGAM once
  # combine DAG_sample rows across the K datasets into one dataset
  combined_sample <- dags[[1]]$DAG_sample
  if (run.K > 1) {
    for (k in 2:run.K) combined_sample <- rbind(combined_sample, dags[[k]]$DAG_sample)
  }
  dags_combined <- list(list(DAG_sample = combined_sample, true_Adjace = true_adjacent[[1]]))
  est_tl <- TLLiNGAM(1, dags_combined, run.hardth, criti.val_c, run.K * run.tuning)$B
  eval_tl <- Evaluation.DAG(est_tl, list(true_adjacent[[1]]), 1)$Eval_result

  ## 4) Single-task TL: run TLLiNGAM separately on each dataset and evaluate across K
  est_single_tl <- vector("list", run.K)
  for (k in seq_len(run.K)) {
    # run TLLiNGAM on dataset k alone
    est_k <- TLLiNGAM(1, list(dags[[k]]), run.hardth, criti.val_c, run.tuning)$B
    # TLLiNGAM(...)$B is a list (length 1 when K=1); extract the matrix if nested
    if (is.list(est_k) && length(est_k) == 1) {
      est_single_tl[[k]] <- est_k[[1]]
    } else {
      est_single_tl[[k]] <- est_k
    }
  }
  eval_single_tl <- Evaluation.DAG(est_single_tl, true_adjacent, run.K)$Eval_result

  ## 5) GGL method (graphical group lasso over the K datasets)
  DAG.sample <- lapply(dags, function(x) x$DAG_sample)
  ggl_results <- JGL(Y = DAG.sample, penalty = "group", lambda1 = 0,
                     lambda2 = 2 * run.tuning / run.K * sqrt(log(max(run.p, run.n)) / run.n),
                     return.whole.theta = TRUE)
  est_ggl <- ggl_results$theta
  # JGL may return theta as a 3D array (p x p x K), a list of matrices, or a single matrix.
  # Normalize into a list of p x p matrices for Evaluation.DAG.
  if (is.array(est_ggl) && length(dim(est_ggl)) == 3) {
    est_ggl_list <- lapply(seq_len(dim(est_ggl)[3]), function(i) as.matrix(est_ggl[,,i]))
  } else if (is.list(est_ggl)) {
    est_ggl_list <- lapply(est_ggl, function(x) as.matrix(x))
  } else if (is.matrix(est_ggl)) {
    est_ggl_list <- list(as.matrix(est_ggl))
  } else {
    stop("Unexpected format for JGL result 'theta' — expected array, list or matrix")
  }
  # Symmetrize each estimated adjacency (make undirected when needed)
  est_ggl_list <- lapply(est_ggl_list, function(mat) { m <- as.matrix(mat); m[lower.tri(m)] <- t(m)[lower.tri(m)]; m })
  eval_ggl <- Evaluation.DAG(est_ggl_list, true_adjacent, run.K)$Eval_result
  
  
  ## 5) put together a small summary table (TPR / FDR / SHD / MCC as in evaluation)
  methods <- c("IntTL", "TL (pooled)", "Single-TL", "GGL")
  results_matrix <- rbind(eval_inttl, eval_tl, eval_single_tl, eval_ggl)
  rownames(results_matrix) <- methods
  return(as.data.frame(results_matrix))
}


#----------------------------
# Run once and show results
#----------------------------
cat("Running a single illustrative simulation...\n")
res <- run_one(run.n, run.p, run.K, run.example, run.pi, run.hardth, run.tuning, criti.val_c)
print(round(res, 3))

