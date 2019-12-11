# Multicore functions / parallel implementations plus speed optimized & utility functions for Seurat 2 & 3 

Multicore and utility functions & implementations for Seurat using doMC / foreach packages.
Implementations are either from me or found on the web. 



## Use case

Some Seurat functions can be fairly slow when run on a single core. To speed up you can use all cores of your computer.

- Seurat 2.x has very limited multicore functionality (ScaleData, Jackstraw). 
- Seurat 3.0 has [implemented multiple functions](https://satijalab.org/seurat/v3.0/future_vignette.html) using _future_.
- Functions here use a `foreach` based parallel implementations/templates are mostly complementary to the implemented to Seurat's implementation 

Tested on OS X, but it is in development.



## Notice

'Future' and 'doMC (foreach)' based parallelisation seem to collide in one case. If you load futures before, using  `NormalizeData` inside a foreach loop, it fails (Error writing to connection).

Solution: do not load  & setup `future` before `NormalizeData`.

```
# After NormalizeData
library(future)
plan("multiprocess", workers = 6)
# So to set Max mem size to 2GB, you would run :
options(future.globals.maxSize = 4000 * 1024^2) 
```



## Content

1. **Seurat.Functions.R**: 
   1. Multi-core / parallelized calculations
   2. Multiplexed plotting / graphics functions without parallelization
2. **Seurat3.Functions.R**
   1. Single-core functions wrapped in multi-core / parallelized foreach loops
3. **Seurat.Multicore.Examples.R**
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
- FilterCells; subset (in v3)
- NormalizeData

### 3. Parallel Implementation by Seurat (3.1)
- NormalizeData
- Jackstraw (from v2)
- ScaleData (from v2)
- FindMarkers
- FindIntegrationAnchors
- FindClusters 

### 4. No Parallel Implementation
- RunTSNE
- RunUMAP

### 5. Other functions implemented / collected here


#### Functions in main script

- `parallel.computing.by.future()`  # Run gc(), load multi-session computing and extend memory limits.
- `seu.Make.Cl.Label.per.cell()`  # Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
- `add.Cl.Label.2.Metadata()`   # Add a metadata column
- `umapNamedClusters()`   # Plot and save umap based on metadata column.
- `clip10Xcellname()`   # Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration).
- `make10Xcellname()`   # Add a suffix
- `read10x()`   # read10x from gzipped and using features.tsv
- `FindAllMarkers.multicore()`  # Multicore version of FindAllMarkers.
- `gene.name.check()`   # Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files.
- `check.genes()`   # Check if genes exist in your dataset.
- `fixZeroIndexing.seurat()`  # Fix zero indexing in seurat clustering, to 1-based indexing
- `getMetadataColumn()` <- mmeta  # Get a metadata column as a named vector
- `getCellIDs.from.meta()`  # Get cellIDs from a metadata column, matching a list of values (using %in%).
- `seu.add.meta.from.table()`   # Add to obj@metadata from an external table
- `seu.PC.var.explained()`  # Determine percent of variation associated with each PC
- `seu.plot.PC.var.explained()`   # Plot the percent of variation associated with each PC

#### Functions in Saving.and.loading.R

- `isave.RDS()` # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
- `isave.RDS.pigz()` # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
- `isave.image()` # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
- `subsetSeuObj.and.Save()` # subset a compressed Seurat Obj and save it in wd.
- `seuSaveRds()` # Save a compressed Seurat Object, with parallel gzip by pgzip
- `sampleNpc()` # Sample N % of a dataframe (obj@metadata), and return the cell IDs.
- `rrRDS()` # Load a list of RDS files with parallel ungzip by pgzip.
- `sssRDS()` #  Save multiple objects into a list of RDS files using parallel gzip by pgzip (optional).
- `ssaveRDS()` # Save an object with parallel gzip by pgzip.
- `rreadRDS()` # Read an object with parallel ungzip by pgzip.
- `snappy_pipe()` # Alternative, fast compression. Low compression rate, lightning fast.
- `pigz_pipe()` # Alternative: normal gzip output (& compression rate), ~*cores faster in zipping.

#### Functions in main script

- multiFeaturePlot.A4
  - multi-core implementation (of generating plots) did not work: it kept hanging at n*100% cpu use.
- multiFeatureHeatmap.A4
- LabelPoint
- LabelUR
- LabelUL
- LabelBR
- LabelBL
- read10x

