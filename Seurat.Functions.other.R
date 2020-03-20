# ------------------------------------------------------------------------------------------
## Seurat.Functions.other.R
# ------------------------------------------------------------------------------------------
# Source: self + web
# source("~/GitHub/Seurat.multicore/Seurat.Functions.other.R")

# ------------------------------


MergeDuplicateGenesSeurat <- function (seu=ls.Seurat[[i]]){ # How many entries are duplicated
  duplicated(rownames(seu))
  if(summarize & y){
    x = table(vec); x= x[x>1]-1;
    print("The following elements have >1 extra copies:")
    print(x) # table formatting requires a separate entry
  }
  return(y)
}


# quick umap ---------------
umap <- function(gene='DLX5', obj =org, pt_size =1) {
  FeaturePlot(object = obj, features.plot = gene, reduction.use = 'umap', pt.size = pt_size)
}

# FeaturePlot with different defaults ------------------------------------------------------------------
aFeaturePlot <- function (object=org, features.plot, min.cutoff = 'q1', max.cutoff = 'q99', 
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



# Save multiple FeaturePlot from a list of genes on A4 jpeg ------------------------
# "Does not work!"
# amultiFeaturePlot.A4 <- function(list.of.genes, object = org, plot.reduction='umap'
#                                  , colors=c("grey", "red"), nr.Col=2, nr.Row =4, cex = ceiling(10/(nr.Col*nr.Row))
#                                  , gene.min.exp = 'q01', gene.max.exp = 'q99'
#                                  , jpeg.res = 225, jpeg.q = 90) {
#   tictoc::tic()
#   list.of.genes = check.genes(list.of.genes, obj = object)
#   lsG = iterBy.over(1:l(list.of.genes), by=nr.Row*nr.Col)
#   
#   ls.plots <- foreach(i=1:l(lsG)) %dopar% {
#     genes = list.of.genes[lsG[[i]]]
#     FeaturePlot(object, features.plot =genes, reduction.use = plot.reduction
#                 , nCol = nr.Col, cols.use = colors, no.axes = T, no.legend = F, vector.friendly = T
#                 , min.cutoff = gene.min.exp, max.cutoff = gene.max.exp, do.return = T
#                 , pt.size = cex)
#   }
#   print ("Plots are calculated.")
#   for (i in 1:l(lsG)) { print(i )
#     plotname = kpp(c(plot.reduction,i, genes, 'jpg' ))
#     jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
#     print(ls.plots[[i]])
#     try.dev.off()
#   }
#   tictoc::toc()
# }; 
# "Does not work!"



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


# ------------------------------------------------------------

if (F) {
  "Very slow for some reason"
  # extended Seurat object
  setClass(
    "xseurat",
    contains="seurat",
    slots=c(p="list", 
            all.genes = "list")
  ) -> xseurat
  # x <- as(org.discarded, "xseurat")
  
  
  extended.Seurat.class <- function(obj) {
    as(obj, "xseurat")
  }
  
  xz = extended.Seurat.class(org)
  is(xz)
  x@p = p
  x@all.genes = all.genes
  is(x)
  
} # if

# ------------------------------------------------------------
# check.genes <- function(list.of.genes = ClassicMarkers, obj = seu3) { # check if genes exist in your dataset
#   missingGenes = setdiff(list.of.genes, rownames(obj))
#   if(length(missingGenes)>0) {iprint("Genes not found in the data:", missingGenes)}
#   intersect(list.of.genes, rownames(obj))
# }



# ------------------------------------------------------------

# ------------------------------------------------------------

# ------------------------------------------------------------

# ------------------------------------------------------------


