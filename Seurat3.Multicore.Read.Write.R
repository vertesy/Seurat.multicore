######################################################################
# Seurat3.Multicore.Read.Write.R
######################################################################
# source ('~/GitHub/Seurat.multicore/Seurat3.Multicore.Read.Write.R')

"Multicore read / write (I/O) functions are https://github.com/vertesy/Seurat.multicore"
"Single core read / write (I/O) functions are in https://github.com/vertesy/Seurat.utils/"

try(irequire(MarkdownReportsDev))
try(irequire(tictoc))

# Save an object -----------------------------------------------
isave.RDS.pigz <- function(object, prefix =NULL, suffix=NULL, showMemObject=T, saveParams =T){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  dir.create(path_rdata)

  if ( "seurat" %in% is(object) & saveParams) {
    try(object@misc$p <- p, silent = T)
    try(object@misc$all.genes  <- all.genes, silent = T)
  }
  if (showMemObject) { memory.biggest.objects() }
  fnameBase = kppu(prefix, substitute(object), suffix, idate())
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",fnameBase , ".Rds")
  tictoc::tic()
  ssaveRDS(object = object, file = fname)
  tictoc::toc()
  say()
}

# seuSaveRds ------------------------------------------------------------------------

seuSaveRds <- function(object = ls.Seurat, tags = setupFlags, use_Original_OutDir = F) { # Save a compressed Seurat Object, with parallel gzip by pgzip
  if (use_Original_OutDir) create_set_Original_OutDir()
  fname.comb.rds = ppp(substitute(object), tags, idate(), ".Rds")
  iprint( fname.comb.rds)
  ssaveRDS(object = object, filename = fname.comb.rds)
  say()
}

# Wrapper layer 2 (top) -----------------------------------------------------------------------------------

# read / load multiple objects
rrRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org"), ...) { # Load a list of RDS files with parallel ungzip by pgzip.
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

# Save multiple objects using pigz by default ---------------------------------------------
sssRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org.ALL", "org"), name.suffix =NULL, ...) { # Save multiple objects into a list of RDS files using parallel gzip by pgzip (optional).
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

### Wrapper layer 1 -----------------------------------------------------------------------------------
ssaveRDS <- function(object, filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], func_type = c("pipe", "builtin")[1], ...) { # Save an object with parallel gzip by pgzip.
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

### rreadRDS -----------------------------------------------------------------------------------
rreadRDS <- function(filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], func_type = c("pipe", "builtin")[1], ...) {  # Read an object with parallel ungzip by pgzip.
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

# snappy_pipe -----------------------------------------------------------------------------------
snappy_pipe <- function(filename, mode="read") { # Alternative, fast compression. Low compression rate, lightning fast.
  if(mode == "read") {
    con <- pipe(paste0("cat ", filename, " | snzip -dc"), "rb")
  } else {
    con <- pipe(paste0("snzip > ", filename), "wb")
  }
  con
}

# pigz_pipe -----------------------------------------------------------------------------------
pigz_pipe <- function(filename, mode="read", cores=6) { # Alternative: normal gzip output (& compression rate), ~*cores faster in zipping.
  if(mode == "read") {
    con <- pipe(paste0("cat ", filename, " | pigz -dcp ", cores), "rb")
  } else {
    con <- pipe(paste0("pigz -p ", cores, " > ", filename), "wb")
  }
  con
}
