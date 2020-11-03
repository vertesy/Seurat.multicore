######################################################################
# Seurat3.Multicore.Generic.Functions.R
######################################################################
# source ('~/GitHub/Packages/Seurat.multicore/Seurat3.Multicore.Generic.Functions.R')

# NOTE:
# Seurat v3 uses the 'future' framework for parallelization.
# https://satijalab.org/seurat/v3.0/future_vignette.html

# These provide an alternative way, advantages/disadvantages are untested.

try(require(Seurat), silent = F)
try(require(doMC), silent = F)
try(require(future), silent = F)

# try(source("~/GitHub/pseudoBulk/barcode.export.from.Seurat/Seurat3.Write.Out.CBCs.for.subset-bam.R"), silent = T)

# ------------------------------------------------------------------------

parallel.computing.by.future <- function(workers_ = 6, maxMemSize = 4000 * 1024^2) { # Run gc(), load multi-session computing and extend memory limits.
  # https://satijalab.org/seurat/v3.0/future_vignette.html
  cat(
    "1. If you load futures before you finished using foreach loops,
    NormalizeData inside a foreach loop fails (Error writing to connection)
    -> I assume 'future' and 'doMC' are not compatible

    2. If you setup computing on e.g. six cores, it runs 6 instances of R with the entire memory space copied.
    If you run out of memory, the system starts using the SSD as memory, and it slows you down extremely extremely extremely.
    -> Therefore it is important to clean up the memory space before setting up multicore computation.

    Loaded: library(future), workers set to 6 (def),set Max mem size to 2GB (def)."   )

  gc(full = T)
  try(memory.biggest.objects())

  library(future)
  # plan("multiprocess", workers = workers_)
  plan("multisession", workers = workers_)
  # So to set Max mem size to 2GB, you would run :
  options(future.globals.maxSize = maxMemSize)
}


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
