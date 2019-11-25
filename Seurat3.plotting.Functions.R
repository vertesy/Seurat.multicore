# ------------------------------------------------------------------------------------------
## Seurat3.plotting.Functions.R
# ------------------------------------------------------------------------------------------
# source("~/GitHub/Seurat.multicore/Seurat3.plotting.Functions.R")

try(require(Seurat), silent = F)
try(require(ggplot2), silent = F)

# ### Functions for plotting
# - sgCellFractionsBarplot


# ------------------------
ssgCellFractionsBarplot.CORE <- function(data, plotname="Cell proportions per genotype") { # sg stands for "seurat ggplot"
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

# ------------------------
sgCellFractionsBarplot.Mseq <- function(data, seedNr=1989, group_by = "genotype",
                                        plotname="Cell proportions") { # sg stands for "seurat ggplot"
  set.seed(seedNr)
  data %>%
    group_by( genotype ) %>% #eval(substitute(group_by))
    sample_n(NrCellsInSmallerDataSet ) %>%
    ssgCellFractionsBarplot.CORE(plotname = plotname)
}

# sgCellFractionsBarplot ------------------------
sgCellFractionsBarplot <- function(data, seedNr=1989, group_by = "orig.ident",
                                   plotname="Cell proportions") { # sg stands for "seurat ggplot"
  set.seed(seedNr)
  data %>%
    group_by( eval(substitute(group_by)) ) %>%
    sample_n(NrCellsInSmallerDataSet ) %>%
    ssgCellFractionsBarplot.CORE(plotname = plotname)
}

# ------------------------
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


# ------------------------
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

# Work in progress ------------------------------------------------------------
