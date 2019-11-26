# ------------------------------------------------------------------------------------------
## Seurat3.plotting.Functions.R
# ------------------------------------------------------------------------------------------
# source("~/GitHub/Seurat.multicore/Seurat3.plotting.Functions.R")

try(require(Seurat), silent = F)
try(require(ggplot2), silent = F)

# ### Functions for plotting
# - qUMAP
# - multiFeaturePlot.A4
# - multiFeatureHeatmap.A4
# - plot.UMAP.tSNE.sidebyside
# - sgCellFractionsBarplot.Mseq
# - ssgCellFractionsBarplot.CORE
# - # sgCellFractionsBarplot
# - sgCellFractionsBarplot
# - ww.variable.exists.and.true
# - # sgCellFractionsBarplot
# - # sgCellFractionsBarplot

# seu.add.parameter.list.2.seurat.object ---------------
qUMAP <- function(f= 'TOP2A', obj =  combined.obj, splitby = NULL, qlow = "q10", qhigh = "q90") { 
  FeaturePlot(combined.obj, reduction = 'umap'
              , min.cutoff = qlow, max.cutoff = qhigh
              , split.by = splitby
              , features = f)
}
# qUMAP(  )


# Save multiple FeaturePlot from a list of genes on A4 jpeg ------------------------
multiFeaturePlot.A4 <- function(list.of.genes, obj = org, plot.reduction='umap', intersectionAssay = c('RNA', 'integrated')[1]
                                , colors=c("grey", "red"), nr.Col=2, nr.Row =4, cex = round(0.1/(nr.Col*nr.Row), digits = 2)
                                , gene.min.exp = 'q01', gene.max.exp = 'q99', subdir =T
                                , jpeg.res = 225, jpeg.q = 90) {
  tictoc::tic()
  ParentDir = OutDir
  if (subdir) create_set_SubDir(... = p0(substitute(list.of.genes),'.', plot.reduction),'/')
  
  list.of.genes = check.genes(list.of.genes = list.of.genes, obj = obj, assay.slot = intersectionAssay)
  lsG = iterBy.over(1:l(list.of.genes), by=nr.Row*nr.Col)
  for (i in 1:l(lsG)) { 
    genes = list.of.genes[lsG[[i]]]
    iprint(i,genes )
    plotname = kpp(c(plot.reduction,i, genes, 'jpg' ))
    
    plot.list = FeaturePlot(object = obj, features =genes, reduction = plot.reduction, combine = F
                            , ncol = nr.Col, cols = colors 
                            , min.cutoff = gene.min.exp, max.cutoff = gene.max.exp
                            , pt.size = cex)
    
    for(i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + NoLegend() + NoAxes()
    }
    
    ggsave(filename = plotname, width = wA4, height = hA4, 
           plot = cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row)
    )
  }
  
  if (subdir) create_set_OutDir(... = ParentDir)
  tictoc::toc()
}; 








# Save multiple FeatureHeatmaps from a list of genes on A4 jpeg -----------------------
# code for quantile: https://github.com/satijalab/seurat/blob/master/R/plotting_internal.R

multiFeatureHeatmap.A4 <- function(list.of.genes, obj = org, gene.per.page=5
                                   , group.cells.by= "batch", plot.reduction='umap'
                                   , cex = iround(3/gene.per.page), sep_scale = F
                                   , gene.min.exp = 'q5', gene.max.exp = 'q95'
                                   , jpeg.res = 225, jpeg.q = 90) {
  
  tictoc::tic()
  list.of.genes = check.genes(list.of.genes, obj = obj)
  
  lsG = iterBy.over(1:l(list.of.genes), by=gene.per.page)
  for (i in 1:l(lsG)) { print(i )
    genes = list.of.genes[lsG[[i]]]
    plotname = kpp(c("FeatureHeatmap",plot.reduction,i, genes, 'jpg' ))
    print(plotname)
    jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
    try(
      FeatureHeatmap(obj, features.plot =genes , group.by = group.cells.by 
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

plot.UMAP.tSNE.sidebyside <- function(obj = org, grouping = 'res.0.6',
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
  
  p1 <- DimPlot(object = obj, reduction.use = "tsne", no.axes = no_axes, cells.use = cells_use
                , no.legend = no_legend, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size) + 
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))
  
  p2 <- DimPlot(object = obj, reduction.use = "umap", no.axes = no_axes, cells.use = cells_use
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


# ------------------------
sgCellFractionsBarplot.Mseq <- function(data, seedNr=1989, group_by = "genotype",
                                        plotname="Cell proportions") { # sg stands for "seurat ggplot"
  set.seed(seedNr)
  data %>%
    group_by( genotype ) %>% #eval(substitute(group_by))
    sample_n(NrCellsInSmallerDataSet ) %>%
    ssgCellFractionsBarplot.CORE(plotname = plotname)
}

# ------------------------
ssgCellFractionsBarplot.CORE <- function(data, plotname="Cell proportions per ...") { # sg stands for "seurat ggplot"
  data %>%
    ggplot( aes(fill=genotype,  #eval(substitute(group_by))
                x = if (ww.variable.exists.and.true(p$'clusternames.are.defined')) Cl.names else integrated_snn_res.0.3)) +
    ggtitle(plotname) +
    geom_bar( position="fill" ) +
    geom_hline(yintercept=.5, color='darkgrey')  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label=..count..), stat='count',position = position_fill(vjust=0.5)) +
    labs(x = "Clusters", y = "Fraction")
}

# # sgCellFractionsBarplot ------------------------
# sgCellFractionsBarplot <- function(data, seedNr=1989, group_by = "orig.ident",
#                                    plotname="Cell proportions") { # sg stands for "seurat ggplot"
#   set.seed(seedNr)
#   data %>%
#     group_by( eval(substitute(group_by)) ) %>%
#     sample_n(NrCellsInSmallerDataSet ) %>%
#     ssgCellFractionsBarplot.CORE(plotname = plotname)
# }

# ------------------------
sgCellFractionsBarplot <- function(data, seedNr=1989, group_by = "orig.ident", fill_by="experiment",
                                   label_sample_count=T, plotname="Cell proportions per ...") { # sg stands for "seurat ggplot"
  # print(fill_by)
  set.seed(seedNr)
  data %>%
    group_by( eval(substitute(group_by)) ) %>%
    sample_n(NrCellsInSmallerDataSet ) %>%
    
    ggplot(aes_string(fill= fill_by)) +
    aes(x=Cl.names) + # ggplot(aes(fill= genotype)) +
    # ggplot(aes(fill= genotype, x=Cl.names)) + #OLD way
    geom_hline(yintercept=.5)  +
    geom_bar( position="fill" ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    { if (label_sample_count) geom_text(aes(label=..count..), stat='count', position = position_fill(vjust=0.5)) } +
    ggtitle(plotname) +
    labs(x = "Clusters", y = "Fraction")
}


# ------------------------
# ------------------------
# ------------------------

# Work in progress ------------------------------------------------------------


#' ww.variable.exists.and.true
#'
#' Check if a variable name is defined, and if so, is it TRUE
#' @param var A variable 
#' @param alt.message Alternative message if the variable + path does not exist. FALSE or string.
#' @export
#' @examples ww.variable.and.path.exists(path = B, alt.message = "Hello, your path/var does not exist.")

ww.variable.exists.and.true <- function(var = al2, alt.message = NULL) {
  Variable.Name = substitute(var)
  if (exists(as.character(Variable.Name))) {
    if (isTRUE(var)) {
      TRUE
    } else {
      cat("Variable", Variable.Name," is not true: ", var)
      FALSE
    }
  } else {
    if (is.null(alt.message) ) {
      iprint("Variable", Variable.Name, "does not exist.")
    } else {
      cat(alt.message)
    }
    FALSE
  }
}
# 
# al1=T; al3=F; al4=3232; # al2 not defined
# ww.variable.exists.and.true(al1)
# ww.variable.exists.and.true(al2)
# ww.variable.exists.and.true(al3)
# ww.variable.exists.and.true(al4)



# ------------------------

save2umaps.A4 <- function(plot_list, pname = F) {
  if (pname ==F) pname = substitute(plot_list)
  p1 = plot_grid(plotlist = plot_list, nrow = 2, ncol = 1, labels = LETTERS[1:l(plot_list)]  )
  save_plot(plot = p1, filename = extPNG(pname), base_height = hA4, base_width = wA4)
}

# ------------------------
save4umaps.A4 <- function(plot_list, pname = F) {
  if (pname==F) pname = substitute(plot_list)
  p1 = plot_grid(plotlist = plot_list, nrow = 2, ncol = 2, labels = LETTERS[1:l(plot_list)]  )
  save_plot(plot = p1, filename = extPNG(pname), base_height = wA4, base_width = hA4)
}


# ------------------------
# ------------------------
# ------------------------
# ------------------------
# 
# sgCellFractionsBarplot <- function(data, seedNr=1989, group_by = "orig.ident") { # sg stands for "seurat ggplot"
#   set.seed(seedNr)
#   data %>%
#     group_by( eval(substitute(group_by)) ) %>%
#     sample_n(NrCellsInSmallerDataSet ) %>%
# 
#     ggplot(aes(fill=experiment, x=Cl.names)) +
#     geom_hline(yintercept=.5)  +
#     geom_bar( position="fill" ) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(x = "Clusters", y = "Fraction")
# }

# sgCellFractionsBarplot <- function(data, seedNr=1989, group_by = "orig.ident", fill_by="experiment",
#                                    label_sample_count=T, plotname="Cell proportions per ...") { # sg stands for "seurat ggplot"
#   # print(fill_by)
#   set.seed(seedNr)
#   data %>%
#     group_by( eval(substitute(group_by)) ) %>%
#     sample_n(NrCellsInSmallerDataSet ) %>%
#     
#     ggplot(aes_string(fill= fill_by)) +
#     aes(x=Cl.names) + # ggplot(aes(fill= genotype)) +
#     # ggplot(aes(fill= genotype, x=Cl.names)) + #OLD way
#     geom_hline(yintercept=.5)  +
#     geom_bar( position="fill" ) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     { if (label_sample_count) geom_text(aes(label=..count..), stat='count', position = position_fill(vjust=0.5)) } +
#     ggtitle(plotname) +
#     labs(x = "Clusters", y = "Fraction")
# }

# Work in progress ------------------------------------------------------------
