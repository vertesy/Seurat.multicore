# ------------------------------------------------------------------------------------------
## Seurat2.Multicore.Functions.R
# ------------------------------------------------------------------------------------------
# Source: self + web

try(require(Seurat), silent = F)
try(require(doMC), silent = F)


# FindAllMarkers.multicore ------------------------------

FindAllMarkers.multicore <- function(obj = org, min_pct = 0.2, logfc_threshold=0.5, only_pos=F, wait=10, resolution='RNA_snn_res.0.5', nCores =6 ){ # Multicore version of FindAllMarkers.
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



# Work in progress ------------------------------------------------------------
