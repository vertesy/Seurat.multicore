# ------------------------------------------------------------------------------------------
## Seurat.Functions.R
# ------------------------------------------------------------------------------------------
# Source: self + web
# for Seurat v2.X

require(Seurat)
require(doMC)
source("~/GitHub/Seurat.multicore/Seurat.Functions.other.R")
source("~/GitHub/Seurat.multicore/Saving.and.loading.R")



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

FindAllMarkers.multicore <- function(obj = org, min_pct = 0.2, logfc_threshold=0.5, only_pos=F, wait=10,resolution='res.0.5', nCores =6 ){
  tictoc::tic()
  nrClusters=length(unique(obj@meta.data[,resolution]))
  N=nrClusters-1

  j=seq(from = 0, by = wait, length.out = nCores)
  j=rep(j,nrClusters)
  ls.DE <- foreach(i=0:N) %dopar% {
    Sys.sleep(j[i+1])
    FindMarkers(obj, ident.1=i, only.pos = only_pos, min.pct=min_pct, logfc.threshold = logfc_threshold)
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
# check.genes ---------------------------------------
check.genes <- function(list.of.genes = ClassicMarkers, obj = org) { # check if genes exist in your dataset
  missingGenes = setdiff(list.of.genes, rownames(obj@data))
  if(length(missingGenes)>0) {iprint("Genes not found in the data:", missingGenes)}
  intersect(list.of.genes, rownames(obj@data))
}


# Save multiple FeaturePlot from a list of genes on A4 jpeg ------------------------
multiFeaturePlot.A4 <- function(list.of.genes, object = org, plot.reduction='umap'
                                , colors=c("grey", "red"), nr.Col=2, nr.Row =4, cex = ceiling(10/(nr.Col*nr.Row))
                                , gene.min.exp = 'q01', gene.max.exp = 'q99'
                                , jpeg.res = 225, jpeg.q = 90) {
  tictoc::tic()
  list.of.genes = check.genes(list.of.genes, obj = object)
  lsG = iterBy.over(1:l(list.of.genes), by=nr.Row*nr.Col)
  for (i in 1:l(lsG)) { print(i )
    genes = list.of.genes[lsG[[i]]]
    plotname = kpp(c(plot.reduction,i, genes, 'jpg' ))

    jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
    # try(
    FeaturePlot(object, features.plot =genes, reduction.use = plot.reduction
                , nCol = nr.Col, cols.use = colors, no.axes = T, no.legend = F, vector.friendly = T
                , min.cutoff = gene.min.exp, max.cutoff = gene.max.exp, do.return = F
                , pt.size = cex)
    # , silent = F    )
    try.dev.off()
  }
  tictoc::toc()
};

# Save multiple FeatureHeatmaps from a list of genes on A4 jpeg -----------------------
# code for quantile: https://github.com/satijalab/seurat/blob/master/R/plotting_internal.R

multiFeatureHeatmap.A4 <- function(list.of.genes, object = org, gene.per.page=5
                                   , group.cells.by= "batch", plot.reduction='umap'
                                   , cex = iround(3/gene.per.page), sep_scale = F
                                   , gene.min.exp = 'q5', gene.max.exp = 'q95'
                                   , jpeg.res = 225, jpeg.q = 90) {

  tictoc::tic()
  list.of.genes = check.genes(list.of.genes, obj = object)

  lsG = iterBy.over(1:l(list.of.genes), by=gene.per.page)
  for (i in 1:l(lsG)) { print(i )
    genes = list.of.genes[lsG[[i]]]
    plotname = kpp(c("FeatureHeatmap",plot.reduction,i, genes, 'jpg' ))
    print(plotname)
    jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
    try(
      FeatureHeatmap(object, features.plot =genes , group.by = group.cells.by
                     , reduction.use = plot.reduction, do.return = F
                     , sep.scale = sep_scale, min.exp = gene.min.exp, max.exp = gene.max.exp
                     , pt.size = cex, key.position = "top")
      , silent = F
    )
    try.dev.off()
  }
  tictoc::toc()
}


# plot.UMAP.tSNE.sidebyside ---------------------------------------------------------------------

plot.UMAP.tSNE.sidebyside <- function(object = org, grouping = 'res.0.6',
                                      no_legend = F,
                                      do_return = TRUE,
                                      do_label = T,
                                      label_size = 10,
                                      vector_friendly = TRUE,
                                      cells_use = NULL,
                                      no_axes = T,
                                      pt_size = 0.5,
                                      name.suffix = NULL,
                                      width = hA4, heigth = 1.75*wA4, filetype = "pdf") { # plot a UMAP and tSNE sidebyside

  p1 <- DimPlot(object = object, reduction.use = "tsne", no.axes = no_axes, cells.use = cells_use
                , no.legend = no_legend, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size) +
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))

  p2 <- DimPlot(object = object, reduction.use = "umap", no.axes = no_axes, cells.use = cells_use
                , no.legend = T, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size) +
    ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

  plots = plot_grid(p1, p2, labels=c("A", "B"), ncol = 2)
  plotname=kpp( 'UMAP.tSNE', grouping, name.suffix, filetype)

  cowplot::save_plot(filename = plotname, plot = plots
                     , ncol = 2 # we're saving a grid plot of 2 columns
                     , nrow = 1 # and 2 rows
                     , base_width = width
                     , base_height = heigth
                     # each individual subplot should have an aspect ratio of 1.3
                     # , base_aspect_ratio = 1.5
  )
}


# replace zero indexed clusternames ------------------------------------------------
fixZeroIndexing.seurat <- function(ColName.metadata = 'res.0.6', obj=org) { # fix zero indexing seurat clustering
  obj@meta.data[ ,ColName.metadata] =  as.numeric(obj@meta.data[ ,ColName.metadata])+1
  print(obj@meta.data[ ,ColName.metadata])
  return(obj)
}


# get Cells from metadata  ------------------------------------------------
mmeta <- function(ColName.metadata = 'batch', obj = org, as_numeric =F) { # get a metadata column as a named vector
  x = as.named.vector(obj@meta.data[ ,ColName.metadata, drop=F])
  if (as_numeric) {
    as.numeric.wNames(x)+1
  } else {x}
}

# GetCellIDs from metadata ---------------
GetCellIDs.from.meta <- function(obj=org, ColName.meta = 'res.0.6', values = 1) {
  idx.matching.cells = which(obj@meta.data[ , ColName.meta] %in% values)
  iprint(l(idx.matching.cells), 'cells found.')
  return(rownames(obj@meta.data)[idx.matching.cells])
}
# GetCellIDs.from.meta()




# Work in progress ------------------------------------------------------------
