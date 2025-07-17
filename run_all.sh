#!/bin/bash
 
echo "Running R scripts..."
Rscript table1_cv.R
Rscript table1-moreEdges.R
Rscript table2_cv.R
Rscript table3_cv.R
Rscript table4_cv.R
Rscript table_unb2.R
Rscript table_unb1.R
 
echo "All tasks completed."
