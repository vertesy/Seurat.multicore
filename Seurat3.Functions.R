# ------------------------------------------------------------------------------------------
## Seurat.Functions.R
# ------------------------------------------------------------------------------------------
# Source: self + web

# NOTE:
# Seurat v3 uses the 'future' framework for parallelization.
# https://satijalab.org/seurat/v3.0/future_vignette.html

# These provide an alternative way, advantages/disadvantages are untested.

try(require(Seurat), silent = F)
try(require(doMC), silent = F)

try(source("~/GitHub/Seurat.multicore/Seurat.Functions.other.R"), silent = T)
try(source("~/GitHub/Seurat.multicore/Saving.and.loading.R"), silent = T)
try(source("~/GitHub/Seurat.multicore/Seurat3.Write.Out.CBCs.for.subset-bam.R"), silent = T)
try(source("~/GitHub/Seurat.multicore/Seurat3.plotting.Functions.R"), silent = T)

# ### Functions
# For parallel processing of other functions see 
# 
# - read10x
# - FindAllMarkers.multicore
# - multiFeaturePlot.A4
# - multiFeatureHeatmap.A4
# - LabelPoint
# - LabelUR
# - LabelUL
# - LabelBR
# - LabelBL



# read10x from gzipped and using features.tsv [from SO]------------------------
read10x <- function(dir) {
  tictoc::tic()
  names <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
  for (i in 1:length(names)) {
    R.utils::gunzip(paste0(dir, "/", names[i], ".gz"))
  }
  file.copy(paste0(dir, "/features.tsv"), paste0(dir, "/genes.tsv"))
  mat <- Seurat::Read10X(dir)
  file.remove(paste0(dir, "/genes.tsv"))
  for (i in 1:length(names)) {
    R.utils::gzip(paste0(dir, "/", names[i]))
  }
  tictoc::toc()
  mat
}

# FindAllMarkers.multicore ------------------------------

FindAllMarkers.multicore <- function(obj = org, min_pct = 0.2, logfc_threshold=0.5, only_pos=F, wait=10, resolution='RNA_snn_res.0.5', nCores =6 ){
  tictoc::tic()
  nrClusters=length(unique(obj@meta.data[,resolution]))
  N=nrClusters-1

  j=seq(from = 0, by = wait, length.out = nCores)
  j=rep(j,nrClusters)
  ls.DE <- foreach(i=0:N) %dopar% {
    Sys.sleep(j[i+1])
    FindMarkers(object = obj, ident.1=i, only.pos = only_pos, min.pct=min_pct, logfc.threshold = logfc_threshold)
  };  
  tictoc::toc()
  return(ls.DE)
}


# FindAllMarkers.multicore2 ------------------------------
# https://github.com/satijalab/seurat/issues/457
# M=length(unique(mydataset@meta.data$ident1stLevel))
# N=M-1
# FindMarker.wrapper <- function(x){
#   FindMarkers(mydataset,ident.1=x, only.pos = TRUE, min.pct=0.25)
# }
# Markers <- bplapply(0:N, FindMarker.wrapper,BPPARAM=MulticoreParam(3))

# ------------------------------

# gene.name.check for read .mtx /write .rds script ---------------------------------------
gene.name.check <- function(Seu.obj = ls.Seurat[[1]] ) {
  rn = rownames(GetAssayData(object = Seu.obj, slot = "counts"))
  llprint("### Gene name pattern")
  
  llogit('`rn = rownames(GetAssayData(object = ls.Seurat[[1]], slot = "counts"))`')
  llogit('`head(grepv(rn, pattern = "-"), 10)`')
  print('pattern = -')
  llprint(head(grepv(rn, pattern = "-"), 10))
  
  llogit('`head(grepv(rn, pattern = "_"), 10)`')
  print('pattern = _')
  llprint(head(grepv(rn, pattern = "_"), 10))
  
  llogit('`head(grepv(rn, pattern = "\\."), 10)`')
  print('pattern = \\.')
  llprint(head(grepv(rn, pattern = "\\."), 10))
  
  llogit('`head(grepv(rn, pattern = "\\.AS[1-9]"), 10)`')
  print('pattern = \\.AS[1-9]')
  llprint(head(grepv(rn, pattern = "\\.AS[1-9]"), 10))
}


# check.genes ---------------------------------------
check.genes <- function(list.of.genes = ClassicMarkers, obj = combined.obj, assay.slot=c('RNA', 'integrated')[1]) { # check if genes exist in your dataset
  all_genes = rownames(GetAssayData(object = obj, assay = assay.slot, slot = "data")); l(all_genes)
  missingGenes = setdiff(list.of.genes, all_genes)
  if(length(missingGenes)>0) {iprint("Genes not found in the data:", missingGenes)}
  intersect(list.of.genes, all_genes)
}


# replace zero indexed clusternames ------------------------------------------------
fixZeroIndexing.seurat <- function(ColName.metadata = 'res.0.6', obj=org) { # fix zero indexing seurat clustering
  obj@meta.data[ ,ColName.metadata] =  as.numeric(obj@meta.data[ ,ColName.metadata])+1  
  print(obj@meta.data[ ,ColName.metadata])
  return(obj)
}


# get Cells from metadata  ------------------------------------------------
getMetadataColumn <- mmeta <- function(ColName.metadata = 'batch', obj = combined.obj, as_numeric =F) { # get a metadata column as a named vector
  stopifnot(ColName.metadata %in% colnames(obj@meta.data))
  
  x = as.named.vector(obj@meta.data[ ,ColName.metadata, drop=F])
  if (as_numeric) {
    as.numeric.wNames(x)+1
  } else {x}
}

# getMetadataColumn <- function(obj=combined.obj, colname = metaD.CL.colname) {
#   stopifnot(colname %in% colnames(combined.obj@meta.data))
#   return(combined.obj[[colname]][,1])
# }

# GetCellIDs from metadata ---------------
getCellIDs.from.meta <- function(obj=org, ColName.meta = 'res.0.6', values = NA) {
  if (is.na(values)) {
    
  }
  
  idx.matching.cells = which(obj@meta.data[ , ColName.meta] %in% values)
  iprint(l(idx.matching.cells), 'cells found.')
  return(rownames(obj@meta.data)[idx.matching.cells])
}
# getCellIDs.from.meta()

# PCA percent of variation associated with each PC ---------------
seu.PC.var.explained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) =1:length(obj@reductions$pca@stdev)
  return(pct)
}

# plot percent of variation associated with each PC ---------------
seu.plot.PC.var.explained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC
  pct <- seu.PC.var.explained(obj)
  wbarplot(pct, ylab= "% of variation explained" , xlab="Principal Components")
  barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex=.5)
}





# Work in progress ------------------------------------------------------------

# check.genes <- function(list.of.genes = ClassicMarkers, obj = seu3) { # check if genes exist in your dataset
#   missingGenes = setdiff(list.of.genes, rownames(obj))
#   if(length(missingGenes)>0) {iprint("Genes not found in the data:", missingGenes)}
#   intersect(list.of.genes, rownames(obj))
# }


