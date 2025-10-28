# IntTL — Integrative learning of linear non-Gaussian DAGs

This repo contains the implementation and scripts used in the
IntTL paper (integrative learning of linear non-Gaussian directed acyclic
graphs).

## What is here
- `main_functions/` — core algorithms and helpers
  - `Integrative_learning.R` — IntTL implementation and internal helpers
  - `Integrative_high.R` — IntTL for high-dimensional graphs
  - `CV_selection.R` — cross-validation utilities
  - `evaluation.R` — evaluation metrics and `Evaluation.DAG()`
- `simulation_example/` — a short, self-contained example you can run in a
  few seconds: `example1.R`, `example2.R` and `DAGs_generate.R`


## Quick start (recommended)
1. Open R and run:

```r
source('simulation_example/example.R')
```

This example simulates one dataset and compares IntTL to several baselines
(pooled TL, single-task TL, group graphical lasso, MD-LiNGAM).
 
## Dependencies
Install required packages in R (minimal list used by examples):

```r
install.packages(c('Matrix', 'igraph', 'JGL', 'gglasso', 'energy'))
```

If a script fails with a missing-package error, install the package and re-run.
 