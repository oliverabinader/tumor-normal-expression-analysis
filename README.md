# tumor-normal-log2fc-manual-calculation
Gene Mean Expression and log2 Fold Change between 2 groups 

This repository contains simple reusable R functions for:

- calculating mean expression values for selected genes across sample groups
- computing log2 fold-change between two groups from log2-normalized expression data
- generating violin plots for tumor vs normal expression comparison

## Notes
The fold-change calculation assumes the expression data is already log2-normalized.  
Therefore, the difference in group means is reported as `log2_fc`.
