# ------------------------------------------------------------------------------------------
## Seurat.Multicore.Examples.R
# ------------------------------------------------------------------------------------------
# DO NOT RUN - these are just code snippets not a ready to go pipeline

# rm(list=ls())
require(tictoc)
require(dplyr)
require(gtools)
require(cowplot)
require(stringr)
require(tibble)
library(Seurat)
source ('~/GitHub/TheCorvinas/R/CodeAndRoll.R')
source(".../Seurat.Functions.R")
require(MarkdownReportsDev)
require(doMC)
registerDoMC(6)

OutDir= "/Users/..."
setup_MarkdownReports(OutDir = OutDir)

p=NULL
p$min.cells = 10
p$min.features = 200
# filtering -----
p$thr.lp.nGene = 6000
p$thr.hp.nGene = 1000
p$thr.lp.mito = 0.15
p$thr.lp.ribo = 0.3


p$'min.pct' = 0.2
p$'logfc.threshold' = 0.5


# Load the dataset
InputDir = ""
Tag = "filtered_feature_bc_matrix."
inNames = c(  "A"
              , "B"
              , "C"
              , "D")

# v.project = str_split_fixed(inNames, '\\.', 3)[,3]
# genotype = str_split_fixed(inNames, '\\.', 4)[,4]
# batch = as.factor(str_split_fixed(v.project, '', 3)[,2])


# read in FILTERED -------------------------------------------------------------------------------------
if (FirstTime) {
  ls.input <- ls.Seurat <- list.fromNames(inNames)
  
  # Read 10x
  ls.input <- foreach(i=1:l(inNames)) %dopar% {  
    read10x(p0(InputDir, Tag, inNames[i])) 
  }; 

  # Create Seurat
  ls.Seurat <- foreach(i=1:l(inNames)) %dopar% {
    CreateSeuratObject(raw.data =  ls.input[[i]], project = v.project[i], 
                       min.cells = p$min.cells, min.features = p$min.features)
  }; 
  
  # Save Seurat
  OutputDir = "/Users/..."
  foreach(i=1:l(inNames)) %dopar% {
    saveRDS(ls.Seurat[[i]], file = ppp(p0(OutputDir, 'filtered.', inNames[i]), 'min.cells', p$min.cells, 'min.features', p$min.features,"Rds") )
  }; sayy()
  
} else {
  # Read Seurat  
  ls.Seurat <- foreach(i=1:l(inNames)) %dopar% {
    readRDS(file = ppp(p0(InputDir, 'filtered.', inNames[i]), 'min.cells.10.min.features.200.Rds') )
  }; 
  
}



# Add metadata ------------------------------

batch = as.factor(1:3)

ls.META<- foreach(i=1:l(ls.Seurat)) %dopar% {
  META = ls.Seurat[[i]]@meta.data
  mito.genes <- grep(pattern = "^MT\\.", x = rownames(x = ls.Seurat[[i]]@data), value = TRUE)
  percent.mito <- Matrix::colSums(ls.Seurat[[i]]@data[mito.genes, ])/Matrix::colSums(ls.Seurat[[i]]@data)
  
  ribo.genes <- grep(pattern = "^RPL|^RPS", x = rownames(x = ls.Seurat[[i]]@data), value = TRUE)
  percent.ribo <- Matrix::colSums(ls.Seurat[[i]]@data[ribo.genes, ])/Matrix::colSums(ls.Seurat[[i]]@data)
  
  META$'percent.mito' <- percent.mito
  META$'percent.ribo' <- percent.ribo
  META$'log10.nUMI' <- as.named.vector(META[,'nUMI', drop=F])
  META$'log10.nGenes' <- as.named.vector(META[,'nGene', drop=F])
  nCells = nrow(META)
  META$'batch'    <- rep(batch[i], nCells)
  META
}

# calculation can take long, and subsettign back only worked in normal for loop AFAIR
for (i in 1:3) { ls.Seurat[[i]]@meta.data <- ls.META[[i]] }


# FilterCells ------------------------------

ls2<- foreach(i=1:l(ls2)) %dopar% {
  FilterCells(object = ls2[[i]], subset.names = c("nGene", "percent.mito", "percent.ribo"),
              low.thresholds = c(p$thr.hp.nGene, -Inf, -Inf),
              high.thresholds = c(p$thr.lp.nGene, p$thr.lp.mito, p$thr.lp.ribo))
}; say()
l(ls2); ls2


# NormalizeData ------------------------------
ls2<- foreach(i=1:l(ls2)) %dopar% {
  NormalizeData(object = ls2[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}; say()
l(ls2); ls20

# FindVariableGenes cannot be parallelized sadly   ------------------------------
tic()
for (i in 1:l(inNames)) { print(i)
  ls2[[i]]  <- FindVariableGenes(object = ls2[[i]], mean.function = ExpMean, dispersion.function = LogVMR,
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5 )
  wplot_save_this(plotname = ppp("Var.genes.",inNames[i]), PNG = T)   
}; say()
toc()

# ScaleData is parallelizedialready  ------------------------------
for (i in 1:l(inNames)) { print(i)
  ls2[[i]]  <- ScaleData(object = ls2[[i]], vars.to.regress = ls.vars2regress[[i]], do.par = T, num.cores = 6)
}; say()


