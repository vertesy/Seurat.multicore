# ------------------------------------------------------------------------------------------
## Seurat.Functions.R
# ------------------------------------------------------------------------------------------
# Source: self + web

require(Seurat)
require(doMC)



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
  mat
}

# FindAllMarkers.multicore ------------------------------

FindAllMarkers.multicore <- function(obj = org, min_pct = 0.2, logfc_threshold=0.5, only_pos=F ){
  nrClusters=length(unique(obj@meta.data$'ident'))
  N=nrClusters-1
  
  ls.DE <- foreach(i=0:N) %dopar% {
    FindMarkers(obj,ident.1=N, only.pos = only_pos, min.pct=min_pct, logfc.threshold = logfc_threshold)
  }; 
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
check.genes <- function(list.of.genes = ClassicMarkers, object = org) { # check if genes exist in your dataset
  missingGenes = setdiff(list.of.genes, rownames(object@data))
  if(length(missingGenes)>0) {iprint("Genes not found in the data:", missingGenes)}
  intersect(list.of.genes, rownames(object@data))
}




# Save multiple FeaturePlot from a list of genes on A4 jpeg ------------------------
multiFeaturePlot.A4 <- function(list.of.genes, obj = org, plot.reduction='umap'
                                , colors=c("grey", "red"), nr.Col=2, nr.Row =4, cex = ceiling(10/(nr.Col*nr.Row))
                                , gene.min.exp = 'q01', gene.max.exp = 'q99'
                                , jpeg.res = 225, jpeg.q = 90) {
  
  list.of.genes = check.genes(list.of.genes, object = obj)
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
}; 




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
  plotname=kpp( 'UMAP.tSNE', grouping, filetype)
  cowplot::save_plot(filename = plotname, plot = plots 
                     , ncol = 2 # we're saving a grid plot of 2 columns
                     , nrow = 1 # and 2 rows
                     , base_width = width
                     , base_height = heigth
                     # each individual subplot should have an aspect ratio of 1.3
                     # , base_aspect_ratio = 1.5
                     )
}

# Save multiple FeatureHeatmaps from a list of genes on A4 jpeg -----------------------
# code for quantile: https://github.com/satijalab/seurat/blob/master/R/plotting_internal.R

multiFeatureHeatmap.A4 <- function(list.of.genes, object = org, gene.per.page=5
                                   , group.cells.by= "batch", plot.reduction='umap'
                                   , cex = iround(3/gene.per.page), sep_scale = F
                                   , gene.min.exp = 'q5', gene.max.exp = 'q95'
                                   , jpeg.res = 225, jpeg.q = 90) {
  
  list.of.genes = check.genes(list.of.genes, object = scObj)
  
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
}

# LabelPoint ------------------------
LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                            label = i, size = text.size)
    plot <- plot + annotate("segment",
                            x = x1 + adj.x.s, 
                            xend = x1, 
                            y = y1 +  adj.y.s, 
                            yend = y1,
                            size = segment.size)
  }
  return(plot)
}

#  ------------------------
LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(
    plot,
    genes,
    exp.mat,
    adj.y.t = adj.u.t,
    adj.x.t = adj.r.t,
    adj.y.s = adj.u.s,
    adj.x.s = adj.r.s,
    ...
  ))
}

#  ------------------------
LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(
    plot,
    genes,
    exp.mat,
    adj.y.t = adj.u.t,
    adj.x.t = -adj.l.t,
    adj.y.s = adj.u.s,
    adj.x.s = -adj.l.s,
    ...
  ))
}

#  ------------------------
LabelBR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(
    plot,
    genes,
    exp.mat,
    adj.y.t = -adj.u.t,
    adj.x.t = adj.l.t,
    adj.y.s = adj.u.s,
    adj.x.s = adj.l.s,
    ...
  ))
}


#  ------------------------
LabelBL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(
    plot,
    genes,
    exp.mat,
    adj.y.t = -adj.u.t,
    adj.x.t = -adj.l.t,
    adj.y.s = adj.u.s,
    adj.x.s = -adj.l.s,
    ...
  ))
}





# FeaturePlot with different defaults ------------------------------------------------------------------
aFeaturePlot <- function (object, features.plot, min.cutoff = 'q1', max.cutoff = 'q99', 
                          dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1
                          , cols.use = c("yellow", "red"), pch.use = 16, overlay = FALSE, do.hover = FALSE, 
                          data.hover = "ident", do.identify = FALSE, reduction.use = "umap", 
                          use.imputed = FALSE, nCol = NULL, no.axes = T, no.legend = F, 
                          coord.fixed = FALSE, dark.theme = FALSE, do.return = FALSE, 
                          vector.friendly = TRUE, png.file = NULL, png.arguments = c(10, 10, 100)) 
{
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  if (is.null(x = nCol)) {
    nCol <- 2
    if (length(x = features.plot) == 1) {
      nCol <- 1
    }
    if (length(x = features.plot) > 6) {
      nCol <- 3
    }
    if (length(x = features.plot) > 9) {
      nCol <- 4
    }
  }
  num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) + 
    1
  if (overlay | do.hover) {
    num.row <- 1
    nCol <- 1
  }
  par(mfrow = c(num.row, nCol))
  dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                              slot = "key")
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- as.data.frame(GetCellEmbeddings(object = object, 
                                               reduction.type = reduction.use, dims.use = c(dim.1, dim.2), 
                                               cells.use = cells.use))
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  data.plot$pt.size <- pt.size
  names(x = data.plot) <- c("x", "y")
  data.use <- t(x = FetchData(object = object, vars.all = features.plot, 
                              cells.use = cells.use, use.imputed = use.imputed))
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    ifelse(test = is.na(x = cutoff), yes = min(data.use[feature, 
                                                        ]), no = cutoff)
  }, cutoff = min.cutoff, feature = features.plot)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    ifelse(test = is.na(x = cutoff), yes = max(data.use[feature, 
                                                        ]), no = cutoff)
  }, cutoff = max.cutoff, feature = features.plot)
  check_lengths = unique(x = vapply(X = list(features.plot, 
                                             min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check_lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  if (overlay) {
    pList <- list(BlendPlot(data.use = data.use, features.plot = features.plot, 
                            data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                            cols.use = cols.use, dim.codes = dim.codes, min.cutoff = min.cutoff, 
                            max.cutoff = max.cutoff, coord.fixed = coord.fixed, 
                            no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme))
  }
  else {
    pList <- mapply(FUN = SingleFeaturePlot, feature = features.plot, 
                    min.cutoff = min.cutoff, max.cutoff = max.cutoff, 
                    coord.fixed = coord.fixed, MoreArgs = list(data.use = data.use, 
                                                               data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                                                               cols.use = cols.use, dim.codes = dim.codes, no.axes = no.axes, 
                                                               no.legend = no.legend, dark.theme = dark.theme, 
                                                               vector.friendly = vector.friendly, png.file = png.file, 
                                                               png.arguments = png.arguments), SIMPLIFY = FALSE)
  }
  if (do.hover) {
    if (length(x = pList) != 1) {
      stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
    }
    if (is.null(x = data.hover)) {
      features.info <- NULL
    }
    else {
      features.info <- FetchData(object = object, vars.all = data.hover)
    }
    return(HoverLocator(plot = pList[[1]], data.plot = data.plot, 
                        features.info = features.info, dark.theme = dark.theme, 
                        title = features.plot))
  }
  else if (do.identify) {
    if (length(x = pList) != 1) {
      stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
    }
    return(FeatureLocator(plot = pList[[1]], data.plot = data.plot, 
                          dark.theme = dark.theme))
  }
  else {
    print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
  }
  ResetPar()
  if (do.return) {
    return(pList)
  }
}



# FeaturePlot version ??? ------------------------------------------------------------------

# 
# FeaturePlot2 <-function (object, features.plot, min.cutoff = NA, max.cutoff = NA, 
#                          dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1, 
#                          cols.use = c("yellow", 
#                                       "red"), pch.use = 16, overlay = FALSE, do.hover = FALSE, 
#                          data.hover = "ident", do.identify = FALSE, reduction.use = "tsne", 
#                          use.imputed = FALSE, nCol = NULL, no.axes = FALSE, no.legend = TRUE, 
#                          coord.fixed = FALSE, dark.theme = FALSE, do.return = FALSE, 
#                          vector.friendly = FALSE, png.file = NULL, png.arguments = c(10, 
#                                                                                      10, 100)) 
# {
#   cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
#   if (is.null(x = nCol)) {
#     nCol <- 2
#     if (length(x = features.plot) == 1) {
#       nCol <- 1
#     }
#     if (length(x = features.plot) > 6) {
#       nCol <- 3
#     }
#     if (length(x = features.plot) > 9) {
#       nCol <- 4
#     }
#   }
#   num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) + 
#     1
#   if (overlay | do.hover) {
#     num.row <- 1
#     nCol <- 1
#   }
#   par(mfrow = c(num.row, nCol))
#   dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
#                               slot = "key")
#   dim.codes <- paste0(dim.code, c(dim.1, dim.2))
#   data.plot <- as.data.frame(GetCellEmbeddings(object = object, 
#                                                reduction.type = reduction.use, dims.use = c(dim.1, dim.2), 
#                                                cells.use = cells.use))
#   x1 <- paste0(dim.code, dim.1)
#   x2 <- paste0(dim.code, dim.2)
#   data.plot$x <- data.plot[, x1]
#   data.plot$y <- data.plot[, x2]
#   data.plot$pt.size <- pt.size
#   names(x = data.plot) <- c("x", "y")
#   data.use <- t(x = FetchData(object = object, vars.all = features.plot, 
#                               cells.use = cells.use, use.imputed = use.imputed))
#   min.cutoff <- mapply(FUN = function(cutoff, feature) {
#     ifelse(test = is.na(x = cutoff), yes = min(data.use[feature, 
#                                                         ]), no = cutoff)
#   }, cutoff = min.cutoff, feature = features.plot)
#   max.cutoff <- mapply(FUN = function(cutoff, feature) {
#     ifelse(test = is.na(x = cutoff), yes = max(data.use[feature, 
#                                                         ]), no = cutoff)
#   }, cutoff = max.cutoff, feature = features.plot)
#   check_lengths = unique(x = vapply(X = list(features.plot, 
#                                              min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
#   if (length(x = check_lengths) != 1) {
#     stop("There must be the same number of minimum and maximum cuttoffs as there are features")
#   }
#   if (overlay) {
#     pList <- list(BlendPlot(data.use = data.use, features.plot = features.plot, 
#                             data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
#                             cols.use = cols.use, dim.codes = dim.codes, min.cutoff = min.cutoff, 
#                             max.cutoff = max.cutoff, coord.fixed = coord.fixed, 
#                             no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme))
#   }
#   else {
#     pList <- mapply(FUN = SingleFeaturePlot, feature = features.plot, 
#                     min.cutoff = min.cutoff, max.cutoff = max.cutoff, 
#                     coord.fixed = coord.fixed, MoreArgs = list(data.use = data.use, 
#                                                                data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
#                                                                cols.use = cols.use, dim.codes = dim.codes, no.axes = no.axes, 
#                                                                no.legend = no.legend, dark.theme = dark.theme, 
#                                                                vector.friendly = vector.friendly, png.file = png.file, 
#                                                                png.arguments = png.arguments), SIMPLIFY = FALSE)
#   }
#   if (do.hover) {
#     if (length(x = pList) != 1) {
#       stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
#     }
#     if (is.null(x = data.hover)) {
#       features.info <- NULL
#     }
#     else {
#       features.info <- FetchData(object = object, vars.all = data.hover)
#     }
#     return(HoverLocator(plot = pList[[1]], data.plot = data.plot, 
#                         features.info = features.info, dark.theme = dark.theme, 
#                         title = features.plot))
#   }
#   else if (do.identify) {
#     if (length(x = pList) != 1) {
#       stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
#     }
#     return(FeatureLocator(plot = pList[[1]], data.plot = data.plot, 
#                           dark.theme = dark.theme))
#   }
#   else {
#     print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
#   }
#   ResetPar()
#   if (do.return) {
#     return(pList)
#   }
# }
# 

# Work in progress ------------------------------------------------------------

if (F) {
  p$'resolutions' = c(0.3, 0.6, 0.9)
  p$'num.ccs' = 15
  
  FindClusters.multicore <- function(obj = org, res = p$'resolutions' ){
    
    ls.Clustering <- foreach(i=0:N) %dopar% {
      FindClusters(obj, reduction.type = "cca.aligned" 
                   , resolution = res , dims.use = 1:p$'num.ccs'
                   , plot.SNN = T, print.output = F)
    }; 
    
    
    # extract data
    
    
    clusterings <- seurat@data.info %>% select(contains("res."))
    head(clusterings)
    return(ls.Clustering)
  }
  
  
}

# replace zero indexed clusternames ------------------------------------------------
fixZeroIndexing.seurat <- function(ColName.metadata = 'res.0.6', obj=org) { # fix zero indexing seurat clustering
  obj@meta.data[ ,ColName.metadata] =  as.numeric(obj@meta.data[ ,ColName.metadata])+1  
  print(obj@meta.data[ ,ColName.metadata])
  return(obj)
}


# get Cells from metadata----
mmeta <- function(ColName.metadata = 'batch', obj = org, as_numeric =F) { # get a metadata column as a named vector
  x = as.named.vector(obj@meta.data[ ,ColName.metadata, drop=F])
  if (as_numeric) {
    as.numeric.wNames(x)+1
  } else {x}
}


# mmeta.add <- function(ColName.metadata = 'res.0.6', newdata.vec=rep(1, nrow(obj@meta.data) ), obj = org, as_numeric =F) {
#   stopifnot(length(newdata.vec) =  nrow(obj@meta.data))
#   obj@meta.data[ ,ColName.metadata] <- newdata.vec
#   }


# Save ------------------------
# requires MarkdownReportsDev (github)
# OutDir="~/"

isave <- function(..., showMemObject=T){ # faster saving of workspace, and compression outside R, 
  path_rdata = paste0("~/Documents/Rdata.files/", basename(OutDir))
  dir.create(path_rdata)
  
  if (showMemObject) { memory.biggest.objects() }
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",idate(),...,".Rdata")
  save.image( file = fname, compress=F)
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  system(paste("gzip", fname),  wait = FALSE) # execute in the background
}

