# Modified from run_lidr_grid_metrics_v2.r (Ben Bright, June 2022)
# to run in LiDAR Dataprep on a single 1-square-kilometer height-normalized
# LAZ point cloud without chunking.

library(lidR)
library(modeest) # for mfv() (mode)
library(moments) # for skewness and kurtosis
library(raster)  # For writeRaster()
library(terra)

# TODO: Make Parameters into CLI args
# TODO: Change from default LZW compression? GDAL options: https://gdal.org/drivers/raster/gtiff.html

# Parameters #

lazdir = "/vagrant/lidr-metrics/data/input"
lazstem = "USGS_LPC_DE_Snds_2013_LAS_2015_14262934020004" # The stem of a classified, height-normalized LAZ file
outdir = "/vagrant/lidr-metrics/data/output"
metrics_func_file = "src/lidr_metrics.r" # Relative path to metrics functions file
cell_res = 10 # Grid metrics cell resolution
align = c(0,0)


# Processing #

# Source functions
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
metrics_func = normalizePath(file.path(script.basename, metrics_func_file))
source(metrics_func)

# Run setup
lazfile <- paste0(lazstem, ".laz")
cat("\n", "Processing: ", lazfile, "\n")
ctg_norm <- readLAScatalog(normalizePath(file.path(lazdir, lazfile), mustWork=FALSE))
opt_filter(ctg_norm) <- "-drop_class 7 9 -drop_z_below -30 -drop_z_above 150"
print(ctg_norm)

# Metrics
print("Generating Metrics")
m <- grid_metrics(ctg_norm, ~Metrics(Z, Intensity, ReturnNumber), res=cell_res, start=align)
names(m) = MetricNames
terra::writeRaster(rast(m), normalizePath(file.path(outdir, paste0(lazstem, "_Metrics.tif")), mustWork=FALSE), gdal=c("INTERLEAVE=BAND"))

# StrataMetrics
print("Generating StrataMetrics")
sm <- grid_metrics(ctg_norm, ~StrataMetrics(Z, ReturnNumber), res=cell_res, start=align)
names(sm) = StrataMetricNames
terra::writeRaster(rast(sm), normalizePath(file.path(outdir, paste0(lazstem, "_StrataMetrics.tif")), mustWork=FALSE), gdal=c("INTERLEAVE=BAND"))

# Canopy metrics
print("Generating canopy metrics")
options <- list(raster_alignment = list(res = cell_res, start = align), chunk_alignment=align)
cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="rumple", .options = options)
writeRaster(cm, normalizePath(file.path(outdir, paste0(lazstem, "_canopy_rumple.tif")), mustWork=FALSE))
cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_mean", .options = options)
writeRaster(cm, normalizePath(file.path(outdir, paste0(lazstem, "_canopy_average_height.tif")), mustWork=FALSE))
cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_sd", .options = options)
writeRaster(cm, normalizePath(file.path(outdir, paste0(lazstem, "_canopy_stddev_height.tif")), mustWork=FALSE))
cm <- catalog_sapply(ctg_norm, CHMmetrics, res=cell_res, start=align, M="chm_max", .options = options)
writeRaster(cm, normalizePath(file.path(outdir, paste0(lazstem, "_canopy_maximum_height.tif")), mustWork=FALSE))

print("Processing complete")
