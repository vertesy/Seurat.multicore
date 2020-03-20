# ------------------------------------------------------------------------------------------
## Seurat3.plotting.Functions.R
# ------------------------------------------------------------------------------------------
# source("~/GitHub/Seurat.multicore/Seurat3.plotting.Functions.R")

try(require(Seurat), silent = F)
try(require(ggplot2), silent = F)

#### Functions in Seurat3.plotting.Functions.R

# - `umapHiLightSel()` # Highlight a set of cells based on clusterIDs provided. 
# - `qUMAP()` # Quick umaps 
# - `multiFeaturePlot.A4()` # Save multiple FeaturePlot from a list of genes on A4 jpeg 
# - `multiFeatureHeatmap.A4()` # Save multiple FeatureHeatmaps from a list of genes on A4 jpeg  
# - `plot.UMAP.tSNE.sidebyside()` # plot a UMAP and tSNE sidebyside 
# - `sgCellFractionsBarplot.Mseq()` # Cell Fractions Barplot for MULTI-seq. sg stands for "seurat ggplot". 
# - `ssgCellFractionsBarplot.CORE()` # Cell Fractions Barplots, basic. sg stands for "seurat ggplot". 
# - `sgCellFractionsBarplot()` # Cell Fractions Barplots. sg stands for "seurat ggplot". 
# - `ww.variable.exists.and.true()` # Check if a variable exists and its value is TRUE. 
# - `save2umaps.A4()` # Save 2 umaps on A4. 
# - `save4umaps.A4()` # Save 4 umaps on A4. 

# ---------------
# ---------------
# ---------------
# ---------------


# ------------------------
# ------------------------
# ------------------------

# Work in progress ------------------------------------------------------------

if (F) {
  "Part of MarkdownreportsDev"
  
  #' ww.variable.exists.and.true
  #'
  #' Check if a variable name is defined, and if so, is it TRUE
  #' @param var A variable 
  #' @param alt.message Alternative message if the variable + path does not exist. FALSE or string.
  #' @export
  #' @examples ww.variable.and.path.exists(path = B, alt.message = "Hello, your path/var does not exist.")
  
  ww.variable.exists.and.true <- function(var = al2, alt.message = NULL) { # Check if a variable exists and its value is TRUE.
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

# # sgCellFractionsBarplot ------------------------
# sgCellFractionsBarplot <- function(data, seedNr=1989, group_by = "orig.ident",
#                                    plotname="Cell proportions") { # sg stands for "seurat ggplot"
#   set.seed(seedNr)
#   data %>%
#     group_by( eval(substitute(group_by)) ) %>%
#     sample_n(NrCellsInSmallerDataSet ) %>%
#     ssgCellFractionsBarplot.CORE(plotname = plotname)
# }


# Work in progress ------------------------------------------------------------

