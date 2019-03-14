# Saving and loading R objects: performance testing
# Modified from: https://rpubs.com/jeffjjohnston/rds_compression 


# Save workspace -----------------------------------------------
# requires MarkdownReportsDev (github) and defining OutDir

isave <- function(..., showMemObject=T){ # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/Rdata.files/", basename(OutDir))
  dir.create(path_rdata)
  
  if (showMemObject) { memory.biggest.objects() }
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",idate(),...,".Rdata")
  save.image( file = fname, compress=F)
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  system(paste("gzip", fname),  wait = FALSE) # execute in the background
}



# Wrapper layer 2 (top) -----------------------------------------------------------------------------------

# read / load multiple objects 
rrRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org"), ...) { # load a list of RDS files with parallel unzip pgzip
  tictoc::tic()
  path_rdata = paste0("~/Documents/Rdata.files/", basename(OutDir))
  dir.exists(path_rdata)
  # lsX <- foreach(obj = list_of_objectnames) %dopar% {
  for (obj in list_of_objectnames) {
    fname = MarkdownReportsDev::kollapse(path_rdata , "/", obj)
    tmp.obj = rreadRDS(filename = fname, ..., cores=2)
    # pigzx multi-core decompression is not so useful, with 2 cores ist actually faster:
    
    assigned.name = strsplit(obj, split = '\\.201*')[[1]][1]
    print(paste(" --> Loaded as object: ",assigned.name))
    assign(x = assigned.name,value = tmp.obj, envir = .GlobalEnv)
  }
  tictoc::toc()
}


# require(clipr)
# clip2clip.vector()
# x = c("org.2019.03.12_10h.Rds", "ls2.2019.03.12_10h.Rds", "ls.Seurat.2019.03.12_10h.Rds")
# rrRDS(x)



# save multiple objects using pigz by default
sssRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org.ALL", "org"), ...) { # parallel save RDS
  tictoc::tic()
  base_name <- character()
  path_rdata = paste0("~/Documents/Rdata.files/", basename(OutDir))
  dir.create(path_rdata)
  for (obj in list_of_objectnames) {
    print(obj)
    if ( "seurat" %in% is(obj)) { 
      obj@misc$p = p 
      obj@misc$all.genes = all.genes
    }
    base_name[i] = paste0(obj, '.', idate(),".Rds", collapse = "")
    fname = paste0(path_rdata , "/", base_name[i], collapse = "")
    ssaveRDS( object = get(obj), filename = fname, ...)
  }
  dput(base_name)
  tictoc::toc()
}

# Wrapper layer 1 -----------------------------------------------------------------------------------
ssaveRDS <- function(object, filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], 
                     func_type = c("pipe", "builtin")[1], ...) {
  if(func_type == "builtin") {
    con <- con_func(filename)
  } else {
    con <- con_func(filename, mode="write")
  }
  on.exit(close(con))
  saveRDS(object, con)
}

rreadRDS <- function(filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], 
                     func_type = c("pipe", "builtin")[1], ...) {
  if(func_type == "builtin") {
    con <- con_func(filename)
  } else {
    con <- con_func(filename, mode="read", ...)
  }
  on.exit(close(con))
  readRDS(con)
}

# Pipes -----------------------------------------------------------------------------------

snappy_pipe <- function(filename, mode="read") { # low compression rate, lightning fast
  if(mode == "read") {
    con <- pipe(paste0("cat ", filename, " | snzip -dc"), "rb")
  } else {
    con <- pipe(paste0("snzip > ", filename), "wb")
  }
  con
}


pigz_pipe <- function(filename, mode="read", cores=6) { # normal gzip compression rate, ~*cores faster in zipping
  if(mode == "read") {
    con <- pipe(paste0("cat ", filename, " | pigz -dcp ", cores), "rb")
  } else {
    con <- pipe(paste0("pigz -p ", cores, " > ", filename), "wb")
  }
  con
}