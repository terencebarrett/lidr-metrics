# Modified from run_lidr_grid_metrics_v2.r (Ben Bright, June 2022)
# to run in LiDAR Dataprep on a single 1-square-kilometer height-normalized
# LAZ point cloud without chunking.

library(lidR)
library(modeest) # for mfv() (mode)
library(moments) # for skewness and kurtosis
library(raster)  # For writeRaster()
library(terra)

# TODO: Parallel process this script -- OR -- just run one file and use Python multiprocess on algorithm
# TODO: Name output acc. to input
# TODO: Change from TIF to CSV output
# TODO: Make Parameters into CLI args


# Parameters #

lazdir = "/vagrant/lidr-metrics/data/input/" # Directory of classified, height-normalized LAZ files
outdir = "/vagrant/lidr-metrics/data/output/"
metrics_func = "/vagrant/lidr-metrics/src/lidr_metrics.r" # Location of metrics functions script
cell_res = 10 # Grid metrics cell resolution
align = c(0,0)


# Processing #

# Setup
print('Setting up run')
source(metrics_func) # source metric functions called by grid_metrics()
ctg_norm <- readLAScatalog(lazdir)
opt_filter(ctg_norm) <- "-drop_class 7 9 -drop_z_below -30 -drop_z_above 150"
print(ctg_norm)

# Metrics
print('Generating Metrics')
m <- grid_metrics(ctg_norm, ~Metrics(Z, Intensity, ReturnNumber), res=cell_res, start=align)
names(m) = MetricNames
terra::writeRaster(rast(m), paste0(outdir,"Metrics.tif"))

# StrataMetrics
print('Generating StrataMetrics')
sm <- grid_metrics(ctg_norm, ~StrataMetrics(Z, ReturnNumber), res=cell_res, start=align)
names(sm) = StrataMetricNames
terra::writeRaster(rast(sm), paste0(outdir,"StrataMetrics.tif"))

# - Canopy metrics
print('Generating canopy metrics')
options <- list(raster_alignment = list(res = cell_res, start = align), chunk_alignment=align)

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="rumple", .options = options)
writeRaster(cm, paste0(outdir, "canopy_rumple.tif"))

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_mean", .options = options)
writeRaster(cm, paste0(outdir, "canopy_average_height.tif"))

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_sd", .options = options)
writeRaster(cm, paste0(outdir, "canopy_stddev_height.tif"))

cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_max", .options = options)
writeRaster(cm, paste0(outdir, "canopy_maximum_height.tif"))
