# Seurat multicore

Multicore and untility functions & implementations for Seurat using doMC / foreach packages.
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
- multiFeatureHeatmap.A4
- LabelPoint
- LabelUR
- LabelUL
- LabelBR
- LabelBL
- read10x

