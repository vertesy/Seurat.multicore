# Saving and loading R objects: performance testing
# Modified from: https://rpubs.com/jeffjjohnston/rds_compression 



# Save an object -----------------------------------------------
isave.RDS <- function(object, prefix =NULL, suffix=NULL, showMemObject=T, saveParams =T){ # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  dir.create(path_rdata)
  
  if (showMemObject) { memory.biggest.objects() }
  if ( "seurat" %in% is(obj) & saveParams) { 
    try(obj@misc$p <- p, silent = T)
    try(obj@misc$all.genes  <- all.genes, silent = T)
  }
  fnameBase = kppu(prefix, substitute(object), suffix, idate())
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",fnameBase , ".Rds")
  tictoc::tic()
  saveRDS(object, file = fname, compress=F)
  tictoc::toc()
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  say()
  system(paste("gzip", fname),  wait = FALSE) # execute in the background
}


# Save an object -----------------------------------------------
isave.RDS.pigz <- function(object, prefix =NULL, suffix=NULL, showMemObject=T, saveParams =T){ # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  dir.create(path_rdata)
  
  if ( "seurat" %in% is(obj) & saveParams) { 
    try(obj@misc$p <- p, silent = T)
    try(obj@misc$all.genes  <- all.genes, silent = T)
  }
  if (showMemObject) { memory.biggest.objects() }
  fnameBase = kppu(prefix, substitute(object), suffix, idate())
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",fnameBase , ".Rds")
  tictoc::tic()
  ssaveRDS(object = object, file = fname)
  tictoc::toc()
  say()
}


# Save workspace -----------------------------------------------
# requires MarkdownReportsDev (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r 

isave.image <- function(..., showMemObject=T, options=c("--force", NULL)[1]){ # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/Rdata.files/", basename(OutDir))
  dir.create(path_rdata)
  
  if (showMemObject) { memory.biggest.objects() }
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",idate(),...,".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  save.image( file = fname, compress=F)
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname),  wait = FALSE) # execute in the background
}

# ------------------------------------------------------------------------

subsetSeuObj.and.Save <- function(obj=ORC, fraction = 0.25 ) {
  cellIDs.keep = sampleNpc(metaDF = obj@meta.data, pc = fraction)
  
  obj_Xpc <- subset(obj, cells = cellIDs.keep) # downsample
  saveRDS(obj_Xpc, compress = TRUE,
          file = ppp(p0(InputDir, 'seu.ORC'), l(cellIDs.keep), 'cells.with.min.features', p$min.features,"Rds" ) )
  say()
}


# ------------------------------------------------------------------------

seuSaveRds <- function(object = ls.Seurat, tags = setupFlags, use_Original_OutDir = F) {
  if (use_Original_OutDir) create_set_Original_OutDir()
  fname.comb.rds = ppp(substitute(object), tags, idate(), ".Rds")
  iprint( fname.comb.rds)
  ssaveRDS(object = object, filename = fname.comb.rds)
  say()
}

# ------------------------------------------------------------------------

sampleNpc <- function(metaDF = MetaData[which(Pass),], pc=0.1) {
  cellIDs = rownames(metaDF)
  nr_cells = floor(l(cellIDs) * pc)
  cellIDs.keep = sample(cellIDs, size = nr_cells, replace = FALSE)
  return(cellIDs.keep)
}

# Wrapper layer 2 (top) -----------------------------------------------------------------------------------

# read / load multiple objects 
rrRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org"), ...) { # load a list of RDS files with parallel unzip pgzip
  tictoc::tic()
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  iprint("Looking for files under:" , path_rdata)
  iprint("(Base on OutDir)")
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
sssRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org.ALL", "org"), name.suffix =NULL, ...) { # parallel save RDS
  tictoc::tic()
  base_name <- character()
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  iprint("Files saved under:" , path_rdata)
  
  try(dir.create(path_rdata))
  for (obj in list_of_objectnames) {
    print(obj)
    if ( "seurat" %in% is(obj)) { 
      obj@misc$p = p 
      obj@misc$all.genes = all.genes
    }
    base_name[i] = paste0(obj, '.', name.suffix, '.', idate(),".Rds", collapse = "")
    fname = paste0(path_rdata , "/", base_name[i], collapse = "")
    ssaveRDS( object = get(obj), filename = fname, ...)
  }
  tictoc::toc()
  dput(base_name)
}

# Wrapper layer 1 -----------------------------------------------------------------------------------
ssaveRDS <- function(object, filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], 
                     func_type = c("pipe", "builtin")[1], ...) {
  tictoc::tic()
  if(func_type == "builtin") {
    con <- con_func(filename)
  } else {
    con <- con_func(filename, mode="write")
  }
  on.exit(close(con))
  saveRDS(object, con)
  tictoc::toc()
}

rreadRDS <- function(filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], 
                     func_type = c("pipe", "builtin")[1], ...) {
  tictoc::tic()
  if(func_type == "builtin") {
    con <- con_func(filename)
  } else {
    con <- con_func(filename, mode="read", ...)
  }
  on.exit(close(con))
  obj=readRDS(con)
  tictoc::toc()
  return(obj)
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