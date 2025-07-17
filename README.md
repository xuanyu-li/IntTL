The program files are for the implementation of our IntTL method proposed in the paper NTEGRATIVE LEARNING OF LINEAR NON-GAUSSIAN DIRECTED ACYCLIC GRAPHS WITH APPLICATION ON MULTI-SOURCE GENE REGULATORY NETWORK ANALYSIS

A complete list of the programs:

CV_selection.R

DAGs_generate.R

evaluation.R

Integrative_high.R

Integrative_learning.R

table1_cv.R

table1-moreEdges.R

table2_cv.R

table3_cv.R

table4_cv.R

table_unb2.R

table_unb1.R

MD-LiNGAM

run_all.sh

The whole simulation results can be reproduced by running
```
$ ./run_all.sh
```

 It will load the functions and R packages, and perform simulation trials of the paper. One can perform simulations with different settings by modifying the parameters. The folder MD-LiNGAM contains the files needed to implement the Mdirect method. To implement this method, you should first install the highDLinGam package by the file highDLingam_1.0.tar.gz in the MD-LiNGAM folder. 
