
library(raster)

source("I:/bcb/projects/cms/lidr_scripts/Metrics.r")



# ---------- Metrics
COR.Metrics = rep(NA,length(MetricNames))
for (i in 1:length(MetricNames)) {

lidr = raster(paste0("D:/output/Metrics/",MetricNames[i],".tif"))
fusion = raster(paste0("D:/Annabella2021/Products/Metrics_30METERS/",MetricNames[i],"_30METERS.asc"))
projection(fusion) = projection(lidr)
fusion = projectRaster(fusion, lidr, method="ngb")

COR.Metrics[i] = cor(lidr[], fusion[], use="complete.obs")
}



# ---------- StrataMetrics
COR.StrataMetrics = rep(NA,length(StrataMetricNames))
for (i in 115:length(StrataMetricNames)) {

lidr = raster(paste0("D:/output/StrataMetrics/",StrataMetricNames[i],".tif"))
fusion = raster(paste0("D:/Annabella2021/Products/StrataMetrics_30METERS/",StrataMetricNames[i],"_30METERS.asc"))
projection(fusion) = projection(lidr)
fusion = projectRaster(fusion, lidr, method="ngb")

COR.StrataMetrics[i] = cor(lidr[], fusion[], use="complete.obs")
}



# ---------- TopoMetrics
COR.TopoMetrics = rep(NA,length(TopoMetricNames))
for (i in 1:length(TopoMetricNames)) {

lidr = raster(paste0("D:/output/TopoMetrics/",TopoMetricNames[i],".tif"))
fusion = raster(paste0("D:/Annabella2021/Products/Metrics_30METERS/",TopoMetricNames[i],"_30METERS.asc"))
projection(fusion) = projection(lidr)
fusion = projectRaster(fusion, lidr, method="ngb")

COR.TopoMetrics[i] = cor(lidr[], fusion[], use="complete.obs")
}


# ---------- CHMMetrics
COR.CHMMetrics = rep(NA,length(CHMMetricNames))
for (i in 1:length(CHMMetricNames)) {

lidr = raster(paste0("D:/output/CanopyMetrics/",CHMMetricNames[i],".tif"))
fusion = raster(paste0("D:/Annabella2021/Products/CanopyMetrics_30METERS/",CHMMetricNames[i],"_30METERS.asc"))
projection(fusion) = projection(lidr)
fusion = projectRaster(fusion, lidr, method="ngb")

COR.CHMMetrics[i] = cor(lidr[], fusion[], use="complete.obs")
}


# So:
# Differences in lower StrataMetrics due to different ground normalization
# Differences in Topo curvature metrics (Jeff Evans and Bob cite the same paper for curvature, but it's different)
# I can't figure out how to do FPV in lidR


