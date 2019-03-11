## Seurat multicore

Multicore and untility functions for Seurat using doMC / foreach packages.
Implementations are either from me or found on the web. Tested on OS X, but in development.



### Content

1. **Seurat.Functions.R**: 
   1. Multi-core / parallelized calculations
   2. Multiplexed plotting / graphics functions without parallelization
2. Seurat.Multicore.Examples.R
   1. Single-core functions  wrapped in multi-core / parallelized loops
3. 