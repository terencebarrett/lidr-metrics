# Script Name:  run_lidr_grid_metrics.r
# Description:  Use lidR package to make lidar grid metrics comparable to those of FUSION
# Author:       Ben Bright
# Date:         Nov 2021
# Note:			https://cran.r-project.org/web/packages/lidR/vignettes/lidR-LAScatalog-engine.html

library(lidR)
library(modeest) # for mfv() (mode)
library(moments) # for skewness and kurtosis

# ------- 1. Edit these variables
lazdir = "D:/LAZ_epsg5070/" # Directory of CLASSIFIED laz files
temp = "D:/temp/" # Directory to store normalized point clouds and metric grid chunks
outdir = "D:/output/" # Directory to store gridmetric rasters
metrics_func = "I:/bcb/projects/cms/lidr_scripts/Metrics.r" # Location of metrics functions script
cores = 3 # Number of cores for LAStools to use
cell_res = 30 # Grid metrics cell resolution
chunk_size = 600 # Chunk size in horizontal units, i.e., 600 means 600 x 600 m if horizontal units are meters, should be multiple of cell_res, chunk size of 1000 was too big for my 16 GB RAM machine
chunk_buffer = 30 # Chunk buffer size
align = c(15,15) # Chunk and grid alignment, default is (0,0), CMS Phase 2 alignment is (15,15)


# -------- 2. Create digital terrain model (this could also be done with lidR, which is likely slower)
# - 2.1 Create 1-m DTM
system(paste0('C:/LAStools/bin/lasindex -i ',lazdir,'*.laz -cores ',cores)) # index laz files
system(paste0('C:/LAStools/bin/blast2dem -i ',lazdir,'*.laz -merged -keep_class 2 8 -o ',outdir,'DTM.tif'))
# - 2.2 Create TopoMetrics
source(metrics_func) # Source TopoMetrics() function
dir.create(paste0(outdir,"TopoMetrics"))
TopoMetrics(paste0(outdir,'DTM.tif'), cell_res, align, paste0(outdir,"/TopoMetrics/")) 


# ------- 3. Normalize to heights above ground (this could also be done with LAStools, which is likely faster)
ctg <- readLAScatalog(lazdir) # create catalog
opt_chunk_size(ctg) <- chunk_size # set chunk size
opt_chunk_buffer(ctg) <- chunk_buffer # set chunk buffer
opt_chunk_alignment(ctg) <- align # set chunk alignment
opt_laz_compression(ctg) <- TRUE # output compressed (laz) files
opt_output_files(ctg) <- paste0(temp,"{ID}_norm") # output laz files and naming convention
ctg_norm <- normalize_height(ctg, tin()) # normalize with TIN on the fly, 2.6 GB of laz took 25 min on my 16 GB RAM machine
system(paste0('C:/LAStools/bin/lasindex -i ',temp,'*.laz -cores ',cores)) # index laz files


# ------- 4. Compute gridmetrics
# - 4.1 Set catalog options for normalized point cloud
# - Filter points classified as 7 (low noise) and 9 (water), and points with normalized heights <-30 or >150 m
# - See readLAS(filter = "-help") for filter options
ctg_norm <- readLAScatalog(temp)
opt_filter(ctg_norm) <- "-drop_class 7 9 -drop_z_below -30 -drop_z_above 150"
opt_chunk_size(ctg_norm) <- chunk_size # set chunk size
opt_chunk_buffer(ctg_norm) <- chunk_buffer # set chunk buffer
opt_chunk_alignment(ctg_norm) <- align # set chunk alignment
source(metrics_func) # source metric functions called by grid_metrics()

# - 4.2 Metrics
opt_output_files(ctg_norm) <- paste0(temp,"{ID}_metrics") # output naming convention
m <- grid_metrics(ctg_norm, ~Metrics(Z, Intensity, ReturnNumber), res=cell_res, start=align)
names(m) = MetricNames
dir.create(paste0(outdir,"Metrics"))
writeRaster(m, filename=paste0(outdir,"/Metrics/",MetricNames), bylayer=TRUE, format="GTiff")
# writeRaster(m, paste0(outdir,"Metrics.tif")) # Above is slow; optionally, output as a multiband raster (band names not preserved for .tif, can use .grd or .envi to preserve band names)

# - 4.3 StrataMetrics
opt_output_files(ctg_norm) <- paste0(temp,"{ID}_strataMetrics") # output naming convention
sm <- grid_metrics(ctg_norm, ~StrataMetrics(Z, ReturnNumber), res=cell_res, start=align)
names(sm) = StrataMetricNames
dir.create(paste0(outdir,"StrataMetrics"))
writeRaster(sm, filename=paste0(outdir,"/StrataMetrics/",StrataMetricNames), bylayer=TRUE, format="GTiff")
# writeRaster(m, paste0(outdir,"StrataMetrics.tif")) # Above is slow; optionally, output as a multiband raster (band names not preserved for .tif, can use .grd or .envi to preserve band names)

# - 4.4 CanopyMetrics
# - CHM
opt_output_files(ctg_norm) <- paste0(temp,"{ID}_chm") # output naming convention
chm <- grid_canopy(ctg_norm, res = 1, pitfree(c(0,2,5,10,15), c(0, 1.5))) # Khosravipour et al. pitfree algorithm
writeRaster(chm, paste0(outdir,"CHM.tif"))

# - CHM metrics
options <- list(raster_alignment = list(res = cell_res, start = align), chunk_alignment=align)
dir.create(paste0(outdir,"CanopyMetrics"))

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="rumple", .options = options)
writeRaster(cm, paste0(outdir,"/CanopyMetrics/","canopy_rumple.tif"))

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_mean", .options = options)
writeRaster(cm, paste0(outdir,"/CanopyMetrics/","canopy_average_height.tif"))

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_sd", .options = options)
writeRaster(cm, paste0(outdir,"/CanopyMetrics/","canopy_stddev_height.tif"))

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_max", .options = options)
writeRaster(cm, paste0(outdir,"/CanopyMetrics/","canopy_maximum_height.tif"))

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_fpv", .options = options)
writeRaster(cm, paste0(outdir,"/CanopyMetrics/","canopy_FPV.tif"))









