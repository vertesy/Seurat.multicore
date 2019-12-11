# ------------------------------------------------------------------------------------------
## Seurat.Functions.R
# ------------------------------------------------------------------------------------------
# Source: self + web

# NOTE:
# Seurat v3 uses the 'future' framework for parallelization.
# https://satijalab.org/seurat/v3.0/future_vignette.html

# These provide an alternative way, advantages/disadvantages are untested.

try(require(Seurat), silent = F)
try(require(doMC), silent = F)

try(source("~/GitHub/Seurat.multicore/Seurat.Functions.other.R"), silent = T)
try(source("~/GitHub/Seurat.multicore/Saving.and.loading.R"), silent = T)
try(source("~/GitHub/pseudoBulk/barcode.export.from.Seurat/Seurat3.Write.Out.CBCs.for.subset-bam.R"), silent = T)
try(source("~/GitHub/Seurat.multicore/Seurat3.plotting.Functions.R"), silent = T)

#### Functions

# - parallel.computing.by.future  # Run gc(), load multi-session computing and extend memory limits.
# - seu.Make.Cl.Label.per.cell  # Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
# - add.Cl.Label.2.Metadata   # Add a metadata column
# - umapNamedClusters   # Plot and save umap based on metadata column.
# - clip10Xcellname   # Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration).
# - make10Xcellname   # Add a suffix
# - read10x   # read10x from gzipped and using features.tsv
# - FindAllMarkers.multicore  # Multicore version of FindAllMarkers.
# - gene.name.check   # Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files.
# - check.genes   # Check if genes exist in your dataset.
# - fixZeroIndexing.seurat  # Fix zero indexing in seurat clustering, to 1-based indexing
# - getMetadataColumn <- mmeta  # Get a metadata column as a named vector
# - getCellIDs.from.meta  # Get cellIDs from a metadata column, matching a list of values (using %in%).
# - seu.add.meta.from.table   # Add to obj@metadata from an external table
# - seu.PC.var.explained  # Determine percent of variation associated with each PC
# - seu.plot.PC.var.explained   # Plot the percent of variation associated with each PC

#### Functions in Saving.and.loading.R

# - `isave.RDS()` # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
# - `isave.RDS.pigz()` # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
# - `isave.image()` # faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
# - `subsetSeuObj.and.Save()` # subset a compressed Seurat Obj and save it in wd.
# - `seuSaveRds()` # Save a compressed Seurat Object, with parallel gzip by pgzip
# - `sampleNpc()` # Sample N % of a dataframe (obj@metadata), and return the cell IDs.
# - `rrRDS()` # Load a list of RDS files with parallel ungzip by pgzip.
# - `sssRDS()` #  Save multiple objects into a list of RDS files using parallel gzip by pgzip (optional).
# - `ssaveRDS()` # Save an object with parallel gzip by pgzip.
# - `rreadRDS()` # Read an object with parallel ungzip by pgzip.
# - `snappy_pipe()` # Alternative, fast compression. Low compression rate, lightning fast.
# - `pigz_pipe()` # Alternative: normal gzip output (& compression rate), ~*cores faster in zipping.

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

# ------------------------------------------------------------------------
seu.Make.Cl.Label.per.cell <- function(TopGenes, clID.per.cell) { # Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
  Cl.names_class= TopGenes[ clID.per.cell ]
  Cl.names_wNr = p0(Cl.names_class,' (',names(Cl.names_class),')')
  return(Cl.names_wNr)
}
# seu.Make.Cl.Label.per.cell(TopGenes = TopGenes.Classic, 
#                            clID.per.cell = getMetadataColumn(ColName.metadata = metaD.CL.colname)  )
# ------------------------------------------------------------------------
add.Cl.Label.2.Metadata <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell=Cl.Label.per.cell ) { # Add a metadata column
  obj@meta.data[, metaD.colname ] = Label.per.cell
  iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:", unique(Label.per.cell))
  return(obj)
}  
# combined.obj <- add.Cl.Label.2.Metadata(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell=Cl.Label.per.cell )

# ------------------------------------------------------------------------

umapNamedClusters <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, ext = "png") { # Plot and save umap based on metadata column.
  fname = ppp("Named.clusters", metaD.colname, ext)
  p.named = 
    DimPlot(obj, reduction = "umap", group.by = metaD.colname, label = T) +
    NoLegend() + 
    ggtitle(metaD.colname)
  save_plot(p.named, filename = fname); p.named
}
# umapNamedClusters(obj = combined.obj, metaD.colname = metaD.colname.labeled)

# ------------------------------------------------------------------------



# ------------------------------------------------------------------------
clip10Xcellname <- function(cellnames) str_split_fixed(cellnames, "_", n = 2)[,1] # Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration).

# ------------------------------------------------------------------------
make10Xcellname <- function(cellnames, suffix="_1") p0(cellnames, suffix) # Add a suffix

# ------------------------------------------------------------------------------------------------
read10x <- function(dir) { # read10x from gzipped and using features.tsv
  tictoc::tic()
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
  tictoc::toc()
  mat
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


# FindAllMarkers.multicore2 ------------------------------
# https://github.com/satijalab/seurat/issues/457
# M=length(unique(mydataset@meta.data$ident1stLevel))
# N=M-1
# FindMarker.wrapper <- function(x){
#   FindMarkers(mydataset,ident.1=x, only.pos = TRUE, min.pct=0.25)
# }
# Markers <- bplapply(0:N, FindMarker.wrapper,BPPARAM=MulticoreParam(3))

# ------------------------------

# gene.name.check for read .mtx /write .rds script ---------------------------------------
gene.name.check <- function(Seu.obj = ls.Seurat[[1]] ) { # Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files.
  rn = rownames(GetAssayData(object = Seu.obj, slot = "counts"))
  llprint("### Gene name pattern")
  
  llogit('`rn = rownames(GetAssayData(object = ls.Seurat[[1]], slot = "counts"))`')
  llogit('`head(grepv(rn, pattern = "-"), 10)`')
  print('pattern = -')
  llprint(head(grepv(rn, pattern = "-"), 10))
  
  llogit('`head(grepv(rn, pattern = "_"), 10)`')
  print('pattern = _')
  llprint(head(grepv(rn, pattern = "_"), 10))
  
  llogit('`head(grepv(rn, pattern = "\\."), 10)`')
  print('pattern = \\.')
  llprint(head(grepv(rn, pattern = "\\."), 10))
  
  llogit('`head(grepv(rn, pattern = "\\.AS[1-9]"), 10)`')
  print('pattern = \\.AS[1-9]')
  llprint(head(grepv(rn, pattern = "\\.AS[1-9]"), 10))
}


# check.genes ---------------------------------------
check.genes <- function(list.of.genes = ClassicMarkers, obj = combined.obj, assay.slot=c('RNA', 'integrated')[1]) { # Check if genes exist in your dataset.
  all_genes = rownames(GetAssayData(object = obj, assay = assay.slot, slot = "data")); l(all_genes)
  missingGenes = setdiff(list.of.genes, all_genes)
  if(length(missingGenes)>0) {iprint("Genes not found in the data:", missingGenes)}
  intersect(list.of.genes, all_genes)
}


# replace zero indexed clusternames ------------------------------------------------
fixZeroIndexing.seurat <- function(ColName.metadata = 'res.0.6', obj=org) { # Fix zero indexing in seurat clustering, to 1-based indexing
  obj@meta.data[ ,ColName.metadata] =  as.numeric(obj@meta.data[ ,ColName.metadata])+1  
  print(obj@meta.data[ ,ColName.metadata])
  return(obj)
}


# get Cells from metadata  ------------------------------------------------
getMetadataColumn <- mmeta <- function(ColName.metadata = 'batch', obj = combined.obj, as_numeric =F) { # Get a metadata column as a named vector
  stopifnot(ColName.metadata %in% colnames(obj@meta.data))
  
  x = as.named.vector(obj@meta.data[ ,ColName.metadata, drop=F])
  if (as_numeric) {
    as.numeric.wNames(x)+1
  } else {x}
}

# GetCellIDs from metadata ---------------
getCellIDs.from.meta <- function(obj=org, ColName.meta = 'res.0.6', values = NA) { # Get cellIDs from a metadata column, matching a list of values (using %in%).
  idx.matching.cells = which(obj@meta.data[ , ColName.meta] %in% values)
  iprint(l(idx.matching.cells), 'cells found.')
  return(rownames(obj@meta.data)[idx.matching.cells])
}
# getCellIDs.from.meta()


# Add to obj@metadata from an external table ------------------------------------------------------------------------
seu.add.meta.from.table <- function(obj = seu.ORC, meta = MetaData.ORC, suffix = ".fromMeta") { # Add to obj@metadata from an external table
  NotFound  = setdiff(colnames(obj), rownames(meta))
  Found     = intersect(colnames(obj), rownames(meta))
  if (l(NotFound)) iprint(l(NotFound), 'cells were not found in meta, e.g.: ', trail(NotFound, N=10))
  
  mCols.new = colnames(meta)
  mCols.old = colnames(obj@meta.data)
  overlap = intersect(mCols.new, mCols.old)
  if (l(overlap)) {
    iprint(l(overlap), 'metadata columns already exist in the seurat object: ', overlap, '. These are tagged as: *', suffix)
    colnames(meta)[overlap] = p0(overlap, suffix)
  }
  mCols.add = colnames(meta)
  obj@meta.data[Found, mCols.add] = meta[ Found,]
  
  return(obj)
} # x=seu.add.meta.from.table()



# PCA percent of variation associated with each PC ---------------
seu.PC.var.explained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) =1:length(obj@reductions$pca@stdev)
  return(pct)
}

# plot percent of variation associated with each PC ---------------
seu.plot.PC.var.explained <- function(obj =  combined.obj) { # Plot the percent of variation associated with each PC
  pct <- seu.PC.var.explained(obj)
  wbarplot(pct, ylab= "% of variation explained" , xlab="Principal Components")
  barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex=.5 )
}



plotTheSoup <- function(CellR.OutputDir = "/Users/abel.vertesy/Dropbox/Abel.IMBA/Data/SoupX_pbmc4k_demo/") {
  library(SoupX)
  # Read In ------------------------
  sc = load10X(CellR.OutputDir, keepDroplets = TRUE)
  # Profiling the soup ------------------------
  sc = estimateSoup(sc)
  
  # Plot top gene's expression ----------------------------------------------------------------
  soupProfile = head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
  soupX.cellfree.RNA.profile = 100 * col2named.vector(soupProfile[,1,drop=F])
  wbarplot(soupX.cellfree.RNA.profile
           , ylab="% Reads in the Soup"
           , sub = paste("Within the", basename(CellR.OutputDir), "dataset")
           , tilted_text = T)
  barplot_label(barplotted_variable = soupX.cellfree.RNA.profile
                , labels = percentage_formatter(soupX.cellfree.RNA.profile/100, digitz = 2)
                , TopOffset = .4, srt = 90, cex=.75)
  
  # Plot summarize expression ----------------------------------------------------------------
  soup.RP.sum   <- colSums(soupProfile[grep('^RPL|^RPS', rownames(soupProfile)),])
  soup.RPL.sum   <- colSums(soupProfile[grep('^RPL', rownames(soupProfile)),])
  soup.RPS.sum   <- colSums(soupProfile[grep('^RPS', rownames(soupProfile)),])
  soup.mito.sum <- colSums(soupProfile[grep('^MT-', rownames(soupProfile)),])
  
  soupProfile.summarized <- rbind(
    'Ribosomal' = soup.RP.sum,
    'Ribosomal.L' = soup.RPL.sum,
    'Ribosomal.S' = soup.RPS.sum,
    'Mitochondial' = soup.mito.sum,
    soupProfile[grep('^RPL|^RPS|^MT-', rownames(soupProfile), invert = T),]
  )
  
  NrColumns2Show  = min(10, nrow(soupProfile.summarized))
  soupX.cellfree.RNA.profile.summarized = 100 * col2named.vector(soupProfile.summarized[1:NrColumns2Show,1,drop=F])
  wbarplot(soupX.cellfree.RNA.profile.summarized
           , ylab="% Reads in the Soup"
           , sub = paste("Within the", basename(CellR.OutputDir), "dataset")
           , tilted_text = T)
  barplot_label(barplotted_variable = soupX.cellfree.RNA.profile.summarized
                , labels = percentage_formatter(soupX.cellfree.RNA.profile.summarized/100, digitz = 2)
                , TopOffset = .5)
  remove("sc")
  detach(SoupX)
} # plotTheSoup

# Work in progress ------------------------------------------------------------

BarplotCellsPerObject <- function(ls.Seu = ls.Seurat, plotname="Nr.Cells.After.Filtering", names=F ) {
  cellCounts = unlapply(ls.Seu, ncol)
  names(cellCounts) = if (l(names) == l(ls.Seurat)) names else names(ls.Seurat) 
  wbarplot(cellCounts, plotname = plotname,tilted_text = T, ylab="Cells")
  barplot_label(cellCounts, TopOffset = 500, w = 4)
}




# Convert10Xfolders ------------------------------
Convert10Xfolders <- function(InputDir, min.cells=10, min.features=200, updateHGNC=T) {
  fin <- list.dirs(InputDir)[-1]
  for (i in 1:l(fin)) { print(fin[i])
    pathIN = fin[i]
    fnameIN = basename(fin[i])
    fnameOUT = ppp(p0(InputDir, 'filtered.', fnameIN), 'min.cells', min.cells, 'min.features', min.features,"Rds")
    x <- Read10X(pathIN)
    seu <- CreateSeuratObject(counts = x, project = fnameIN,
                              min.cells = min.cells, min.features = min.features)
    # update----
    if (updateHGNC) seu <- UpdateGenesSeurat(seu)
    saveRDS(seu, file = fnameOUT)
  }
}
# Convert10Xfolders(InputDir = InputDir)


# LoadAllSeurats -------
LoadAllSeurats <- function(InputDir) {
  fin <- list.files(InputDir, include.dirs = F, pattern = "*.Rds")
  ls.Seu <- list.fromNames(fin)
  for (i in 1:l(fin)) {print(fin[i]); ls.Seu[[i]] <- rreadRDS(p0(InputDir, fin[i]))}
  return(ls.Seu)
}
# ls.Seu <- LoadAllSeurats(InputDir = InputDir)



# updateHGNC helper -------
RenameGenesSeurat <- function(SeuObj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]) {
  print("Run this before integration. It only changes SeuObj@assays$RNA@counts, @data and @scale.data")
  RNA <- SeuObj@assays$RNA
  
  if (nrow(RNA) == nrow(newnames)) {
    if(l(RNA@counts)) RNA@counts@Dimnames[[1]] <-         newnames$Suggested.Symbol
    if(l(RNA@data)) RNA@data@Dimnames[[1]] <-             newnames$Suggested.Symbol
    if(l(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames$Suggested.Symbol
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  SeuObj@assays$RNA <- RNA
  return(SeuObj)
}

# updateHGNC -------
UpdateGenesSeurat <- function(seu, species_="human") {
  HGNC.updated <- checkGeneSymbols(rownames(seu), unmapped.as.na = FALSE, map = NULL, species = species_)
  seu <- RenameGenesSeurat(seu, newnames = HGNC.updated)
}


# updateHGNC plot -------
plot.UpdateStats <- function(genes = HGNC.updated[[i]]) {
  (MarkedAsUpdated <- genes[genes$Approved == FALSE, ])
  (AcutallyUpdated <- sum(MarkedAsUpdated[,1] != MarkedAsUpdated[,3]))
  (UpdateStats = c((AcutallyUpdated / nrow(genes)), AcutallyUpdated, nrow(genes)))
  return(UpdateStats)
}



# Work in progress ------------------------------------------------------------

# check.genes <- function(list.of.genes = ClassicMarkers, obj = seu3) { # check if genes exist in your dataset
#   missingGenes = setdiff(list.of.genes, rownames(obj))
#   if(length(missingGenes)>0) {iprint("Genes not found in the data:", missingGenes)}
#   intersect(list.of.genes, rownames(obj))
# }



# getMetadataColumn <- function(obj=combined.obj, colname = metaD.CL.colname) {
#   stopifnot(colname %in% colnames(combined.obj@meta.data))
#   return(combined.obj[[colname]][,1])
# }