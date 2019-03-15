# Seurat multicore

Multicore and utility functions & implementations for Seurat using doMC / foreach packages.
Implementations are either from me or found on the web. 

Tested on OS X, but it is in development.



## Content

1. **Seurat.Functions.R**: 
   1. Multi-core / parallelized calculations
   2. Multiplexed plotting / graphics functions without parallelization
2. **Seurat.Multicore.Examples.R**
   1. Single-core functions wrapped in multi-core / parallelized foreach loops




## Implementations

### 1.Parallel Implementation
- FindAllMarkers.multicore
- rrRDS: read a list of objects
  - rreadRDS: optionally parallel decompression of saved object by pigz_pipe
- sssRDS: save a list of objects
  - ssaveRDS: parallel compression of saved objects by pigz_pipe
    - snappy_pipe: fast single-core, loose compression
    - pigz_pipe

### 2.Parallel templates with foreach
- CreateSeuratObject
- saveRDS
- readRDS
- METADATA
- FilterCells
- NormalizeData

### 3. Parallel Implementation by Seurat
- ScaleData
- Jackstraw

### 4. No Parallel Implementation
- RunTSNE
- RunUMAP

### 5. Other functions implemented / collected here
- multiFeaturePlot.A4
  - multi-core implementation (of generating plots) did not work: it kept hanging at n*100% cpu use.
- multiFeatureHeatmap.A4
- LabelPoint
- LabelUR
- LabelUL
- LabelBR
- LabelBL
- read10x

