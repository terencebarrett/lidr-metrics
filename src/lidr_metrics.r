# Modified from:
#
# Script Name:  Metrics_v2.r
# Description:  Metrics functions to call with lidR
# Author:       Ben Bright
# Date:         Nov 2021, updated June 2022
# Note:			Added some metrics based on input from Terry Barrett at UVM and Bob McGaughey (21 June 2022 email from Andy)

# ------- TopoMetrics based on digital elevation model
TopoMetrics = function(dem_in, res, align, outdir) {

	library(raster)
	library(spatialEco)

	dem = raster(dem_in)
	XMIN = xmin(dem)-xmin(dem)%%res-res-align[1] # Calculate new XMIN so aggregated grid aligns with other metric grids
	YMAX = ymax(dem)-ymax(dem)%%res+res+align[2] # Calculate new YMAX so aggregated grid aligns with other metric grids
	e = extent(XMIN, xmax(dem), ymin(dem), YMAX)
	dem = extend(dem, e)
	dem = aggregate(dem, res) # Aggregate 1-m DEM to resolution

	# --- Get latitude
	latMin = ymin(extent(projectExtent(dem, "+proj=longlat")))
	latMax = ymax(extent(projectExtent(dem, "+proj=longlat")))
	latitude = mean(c(latMin, latMax))

	# --- Elevation
	writeRaster(dem, paste0(outdir,"TOPO_elevation.tif"))

	# --- Aspect
	aspect = terrain(dem, opt='aspect', unit='degrees', neighbors=8)
	writeRaster(aspect, paste0(outdir,"TOPO_aspect.tif"))

	# --- Slope
	slope = terrain(dem, opt='slope', unit='degrees', neighbors=8)
	writeRaster(slope, paste0(outdir,"TOPO_slope.tif"))

	# --- Curvature
	# From Zevenbergen, L.W. & C.R. Thorne (1987). Quantitative Analysis of Land Surface Topography.
	# Earth Surface Processes and Landforms, 12:47-56.
	# Implemented by Evans in spatialEco package
	curv.planform = curvature(dem, type="planform")
	curv.profile = curvature(dem, type="profile")
	curv.total = curvature(dem, type="total")
	writeRaster(curv.planform*100, paste0(outdir, "TOPO_plancurv.tif"))	# FUSION mutiplies curvature by 100,
	writeRaster(curv.profile*100, paste0(outdir, "TOPO_profilecurv.tif"))	# but still big diffs between this and FUSION
	writeRaster(curv.total*100, paste0(outdir, "TOPO_curvature.tif"))

	# --- Solar radiation index
	# From Keating et al. 2007 (see p. 87 of FUSION manual)
	sri = 1 + cos(latitude*(pi/180)) * cos(slope*(pi/180)) + sin(latitude*(pi/180)) * sin(slope*(pi/180)) * cos((180-aspect)*(pi/180))
	writeRaster(sri, paste0(outdir, "TOPO_sri.tif"))

}

TopoMetricNames = c("TOPO_aspect","TOPO_curvature","TOPO_elevation","TOPO_plancurv","TOPO_profilecurv","TOPO_slope","TOPO_sri")

# ------- Metrics (>2 m aboveground)
Metrics = function(z, i, r) {

	twoPlus = which(z > 2)  # all returns >2 m above ground
	first = which(z > 2 & r == 1) # first returns >2 m above ground

  metrics = list(
	ALL_RETURNS_cnt = length(z), # ADDED JUNE 2022
	ALL_RETURNS_all_cnt_2plus = length(z[twoPlus]),
	ALL_RETURNS_elev_AAD_2plus = mean( abs(z[twoPlus] - mean(z[twoPlus])) ),
	ALL_RETURNS_elev_ave_2plus = mean(z[twoPlus]),
	ALL_RETURNS_elev_canopy_relief_ratio = (mean(z[twoPlus]) - min(z[twoPlus])) / (max(z[twoPlus]) - min(z[twoPlus])),
	ALL_RETURNS_elev_CV_2plus = sd(z[twoPlus])/mean(z[twoPlus]),
	ALL_RETURNS_elev_IQ_2plus = IQR(z[twoPlus]),
	ALL_RETURNS_elev_kurtosis_2plus = kurtosis(z[twoPlus]),
	ALL_RETURNS_elev_max_2plus = max(z[twoPlus]),
	ALL_RETURNS_elev_min_2plus = min(z[twoPlus]), # ADDED JUNE 2022
	ALL_RETURNS_elev_mode_2plus = as.numeric(mfv(z[twoPlus])[1]), # ADDED JUNE 2022
	ALL_RETURNS_elev_P01_2plus = quantile(z[twoPlus], 0.01),
	ALL_RETURNS_elev_P05_2plus = quantile(z[twoPlus], 0.05),
	ALL_RETURNS_elev_P10_2plus = quantile(z[twoPlus], 0.10), # ADDED JUNE 2022
	ALL_RETURNS_elev_P20_2plus = quantile(z[twoPlus], 0.20), # ADDED JUNE 2022
	ALL_RETURNS_elev_P25_2plus = quantile(z[twoPlus], 0.25),
	ALL_RETURNS_elev_P30_2plus = quantile(z[twoPlus], 0.30), # ADDED JUNE 2022
	ALL_RETURNS_elev_P40_2plus = quantile(z[twoPlus], 0.40), # ADDED JUNE 2022
	ALL_RETURNS_elev_P50_2plus = quantile(z[twoPlus], 0.50),
	ALL_RETURNS_elev_P60_2plus = quantile(z[twoPlus], 0.60), # ADDED JUNE 2022
	ALL_RETURNS_elev_P70_2plus = quantile(z[twoPlus], 0.70), # ADDED JUNE 2022
	ALL_RETURNS_elev_P75_2plus = quantile(z[twoPlus], 0.75),
	ALL_RETURNS_elev_P80_2plus = quantile(z[twoPlus], 0.80), # ADDED JUNE 2022
	ALL_RETURNS_elev_P90_2plus = quantile(z[twoPlus], 0.90), # ADDED JUNE 2022
	ALL_RETURNS_elev_P95_2plus = quantile(z[twoPlus], 0.95),
	ALL_RETURNS_elev_P99_2plus = quantile(z[twoPlus], 0.99),
	ALL_RETURNS_elev_skewness_2plus = skewness(z[twoPlus]),
	ALL_RETURNS_elev_stddev_2plus = sd(z[twoPlus]),
	ALL_RETURNS_elev_variance_2plus = var(z[twoPlus]),
	ALL_RETURNS_int_AAD_2plus = mean( abs(i[twoPlus] - mean(i[twoPlus])) ),
	ALL_RETURNS_int_ave_2plus = mean(i[twoPlus]),
	ALL_RETURNS_int_CV_2plus = sd(i[twoPlus])/mean(i[twoPlus]),
	ALL_RETURNS_int_IQ_2plus = IQR(i[twoPlus]),
	ALL_RETURNS_int_kurtosis_2plus = kurtosis(i[twoPlus]),
	ALL_RETURNS_int_max_2plus = as.numeric(max(i[twoPlus])),
	ALL_RETURNS_int_min_2plus = as.numeric(min(i[twoPlus])),
	ALL_RETURNS_int_mode_2plus = as.numeric(mfv(i[twoPlus])[1]), # Added by TCB July 2022
	ALL_RETURNS_int_P01_2plus = quantile(i[twoPlus], 0.01),
	ALL_RETURNS_int_P05_2plus = quantile(i[twoPlus], 0.05),
	ALL_RETURNS_int_P10_2plus = quantile(i[twoPlus], 0.10), # ADDED JUNE 2022
	ALL_RETURNS_int_P20_2plus = quantile(i[twoPlus], 0.20), # ADDED JUNE 2022
	ALL_RETURNS_int_P25_2plus = quantile(i[twoPlus], 0.25),
	ALL_RETURNS_int_P30_2plus = quantile(i[twoPlus], 0.30), # ADDED JUNE 2022
	ALL_RETURNS_int_P40_2plus = quantile(i[twoPlus], 0.40), # ADDED JUNE 2022
	ALL_RETURNS_int_P50_2plus = quantile(i[twoPlus], 0.50),
	ALL_RETURNS_int_P60_2plus = quantile(i[twoPlus], 0.60), # ADDED JUNE 2022
	ALL_RETURNS_int_P70_2plus = quantile(i[twoPlus], 0.70), # ADDED JUNE 2022
	ALL_RETURNS_int_P75_2plus = quantile(i[twoPlus], 0.75),
	ALL_RETURNS_int_P80_2plus = quantile(i[twoPlus], 0.80), # ADDED JUNE 2022
	ALL_RETURNS_int_P90_2plus = quantile(i[twoPlus], 0.90), # ADDED JUNE 2022
	ALL_RETURNS_int_P95_2plus = quantile(i[twoPlus], 0.95),
	ALL_RETURNS_int_P99_2plus = quantile(i[twoPlus], 0.99),
	ALL_RETURNS_int_skewness_2plus = skewness(i[twoPlus]),
	ALL_RETURNS_int_stddev_2plus = sd(i[twoPlus]),
	ALL_RETURNS_int_variance_2plus = var(i[twoPlus]),
	FIRST_RETURNS_cnt = length(z[which(r == 1)]), # ADDED JUNE 2022
	FIRST_RETURNS_all_cnt_2plus = length(z[first]),
	FIRST_RETURNS_elev_AAD_2plus = mean( abs(z[first] - mean(z[first])) ),
	FIRST_RETURNS_elev_ave_2plus = mean(z[first]),
	FIRST_RETURNS_elev_canopy_relief_ratio = (mean(z[first]) - min(z[first])) / (max(z[first]) - min(z[first])),
	FIRST_RETURNS_elev_CV_2plus = sd(z[first])/mean(z[first]),
	FIRST_RETURNS_elev_IQ_2plus = IQR(z[first]),
	FIRST_RETURNS_elev_kurtosis_2plus = kurtosis(z[first]),
	FIRST_RETURNS_elev_max_2plus = max(z[first]),
	FIRST_RETURNS_elev_min_2plus = min(z[first]), # ADDED JUNE 2022
	FIRST_RETURNS_elev_mode_2plus = as.numeric(mfv(z[first])[1]), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P01_2plus = quantile(z[first], 0.01),
	FIRST_RETURNS_elev_P05_2plus = quantile(z[first], 0.05),
	FIRST_RETURNS_elev_P10_2plus = quantile(z[first], 0.10), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P20_2plus = quantile(z[first], 0.20), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P25_2plus = quantile(z[first], 0.25),
	FIRST_RETURNS_elev_P30_2plus = quantile(z[first], 0.30), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P40_2plus = quantile(z[first], 0.40), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P50_2plus = quantile(z[first], 0.50),
	FIRST_RETURNS_elev_P60_2plus = quantile(z[first], 0.60), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P70_2plus = quantile(z[first], 0.70), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P75_2plus = quantile(z[first], 0.75),
	FIRST_RETURNS_elev_P80_2plus = quantile(z[first], 0.80), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P90_2plus = quantile(z[first], 0.90), # ADDED JUNE 2022
	FIRST_RETURNS_elev_P95_2plus = quantile(z[first], 0.95),
	FIRST_RETURNS_elev_P99_2plus = quantile(z[first], 0.99),
	FIRST_RETURNS_elev_skewness_2plus = skewness(z[first]),
	FIRST_RETURNS_elev_stddev_2plus = sd(z[first]),
	FIRST_RETURNS_elev_variance_2plus = var(z[first]),
	FIRST_RETURNS_int_AAD_2plus = mean( abs(i[first] - mean(i[first])) ),
	FIRST_RETURNS_int_ave_2plus = mean(i[first]),
	FIRST_RETURNS_int_CV_2plus = sd(i[first])/mean(i[first]),
	FIRST_RETURNS_int_IQ_2plus = IQR(i[first]),
	FIRST_RETURNS_int_kurtosis_2plus = kurtosis(i[first]),
	FIRST_RETURNS_int_max_2plus = as.numeric(max(i[first])),
	FIRST_RETURNS_int_min_2plus = as.numeric(min(i[first])),
	FIRST_RETURNS_int_mode_2plus = as.numeric(mfv(i[first])[1]),
	FIRST_RETURNS_int_P01_2plus = quantile(i[first], 0.01),
	FIRST_RETURNS_int_P05_2plus = quantile(i[first], 0.05),
	FIRST_RETURNS_int_P10_2plus = quantile(i[first], 0.10), # ADDED JUNE 2022
	FIRST_RETURNS_int_P20_2plus = quantile(i[first], 0.20), # ADDED JUNE 2022
	FIRST_RETURNS_int_P25_2plus = quantile(i[first], 0.25),
	FIRST_RETURNS_int_P30_2plus = quantile(i[first], 0.30), # ADDED JUNE 2022
	FIRST_RETURNS_int_P40_2plus = quantile(i[first], 0.40), # ADDED JUNE 2022
	FIRST_RETURNS_int_P50_2plus = quantile(i[first], 0.50),
	FIRST_RETURNS_int_P60_2plus = quantile(i[first], 0.60), # ADDED JUNE 2022
	FIRST_RETURNS_int_P70_2plus = quantile(i[first], 0.70), # ADDED JUNE 2022
	FIRST_RETURNS_int_P75_2plus = quantile(i[first], 0.75),
	FIRST_RETURNS_int_P80_2plus = quantile(i[first], 0.80), # ADDED JUNE 2022
	FIRST_RETURNS_int_P90_2plus = quantile(i[first], 0.90), # ADDED JUNE 2022
	FIRST_RETURNS_int_P95_2plus = quantile(i[first], 0.95),
	FIRST_RETURNS_int_P99_2plus = quantile(i[first], 0.99),
	FIRST_RETURNS_int_skewness_2plus = skewness(i[first]),
	FIRST_RETURNS_int_stddev_2plus = sd(i[first]),
	FIRST_RETURNS_int_variance_2plus = var(i[first]),
	Pct1stRtns_above_2 = length(z[first]) / length(z[which(r == 1)]) * 100,
	Pct1stRtns_above_mean = length(z[which(r == 1 & z > mean(z[first]))]) / length(z[which(r == 1)]) * 100,
	Pct1stRtns_above_mode = length(z[which(r == 1 & z > as.numeric(mfv(z[first])[1]))]) / length(z[which(r == 1)]) * 100, # ADDED JUNE 2022
	PctAllRtns_above_2 = length(z[twoPlus]) / length(z) * 100,
	PctAllRtns_above_mean = length(z[which(z > mean(z[twoPlus]))]) / length(z) * 100,
	PctAllRtns_above_mode = length(z[which(z > as.numeric(mfv(z[twoPlus])[1]))]) / length(z) * 100, # ADDED JUNE 2022
	AllRtns_above_2_dividedby_Tot1stRtns_times100 = length(z[twoPlus]) / length(z[which(r == 1)]) * 100,  # Added by TCB July 2022
	AllRtns_above_mean_dividedby_Tot1stRtns_times100 = length(z[which(z > mean(z[twoPlus]))]) / length(z[which(r == 1)]) * 100,  # Added by TCB July 2022
	AllRtns_above_mode_dividedby_Tot1stRtns_times100 = length(z[which(z > as.numeric(mfv(z[twoPlus])[1]))]) / length(z[which(r == 1)]) * 100  # Added by TCB July 2022
	)

   return(metrics)
}

MetricNames = c(
	"ALL_RETURNS_cnt", # ADDED JUNE 2022
	"ALL_RETURNS_all_cnt_2plus",
	"ALL_RETURNS_elev_AAD_2plus",
	"ALL_RETURNS_elev_ave_2plus",
	"ALL_RETURNS_elev_canopy_relief_ratio",
	"ALL_RETURNS_elev_CV_2plus",
	"ALL_RETURNS_elev_IQ_2plus",
	"ALL_RETURNS_elev_kurtosis_2plus",
	"ALL_RETURNS_elev_max_2plus",
	"ALL_RETURNS_elev_min_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_mode_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P01_2plus",
	"ALL_RETURNS_elev_P05_2plus",
	"ALL_RETURNS_elev_P10_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P20_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P25_2plus",
	"ALL_RETURNS_elev_P30_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P40_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P50_2plus",
	"ALL_RETURNS_elev_P60_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P70_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P75_2plus",
	"ALL_RETURNS_elev_P80_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P90_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_elev_P95_2plus",
	"ALL_RETURNS_elev_P99_2plus",
	"ALL_RETURNS_elev_skewness_2plus",
	"ALL_RETURNS_elev_stddev_2plus",
	"ALL_RETURNS_elev_variance_2plus",
	"ALL_RETURNS_int_AAD_2plus",
	"ALL_RETURNS_int_ave_2plus",
	"ALL_RETURNS_int_CV_2plus",
	"ALL_RETURNS_int_IQ_2plus",
	"ALL_RETURNS_int_kurtosis_2plus",
	"ALL_RETURNS_int_max_2plus",
	"ALL_RETURNS_int_min_2plus",
	"ALL_RETURNS_int_mode_2plus", # Added by TCB July 2022
	"ALL_RETURNS_int_P01_2plus",
	"ALL_RETURNS_int_P05_2plus",
	"ALL_RETURNS_int_P10_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_int_P20_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_int_P25_2plus",
	"ALL_RETURNS_int_P30_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_int_P40_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_int_P50_2plus",
	"ALL_RETURNS_int_P60_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_int_P70_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_int_P75_2plus",
	"ALL_RETURNS_int_P80_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_int_P90_2plus", # ADDED JUNE 2022
	"ALL_RETURNS_int_P95_2plus",
	"ALL_RETURNS_int_P99_2plus",
	"ALL_RETURNS_int_skewness_2plus",
	"ALL_RETURNS_int_stddev_2plus",
	"ALL_RETURNS_int_variance_2plus",
	"FIRST_RETURNS_cnt", # ADDED JUNE 2022
	"FIRST_RETURNS_all_cnt_2plus",
	"FIRST_RETURNS_elev_AAD_2plus",
	"FIRST_RETURNS_elev_ave_2plus",
	"FIRST_RETURNS_elev_canopy_relief_ratio",
	"FIRST_RETURNS_elev_CV_2plus",
	"FIRST_RETURNS_elev_IQ_2plus",
	"FIRST_RETURNS_elev_kurtosis_2plus",
	"FIRST_RETURNS_elev_max_2plus",
	"FIRST_RETURNS_elev_min_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_mode_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P01_2plus",
	"FIRST_RETURNS_elev_P05_2plus",
	"FIRST_RETURNS_elev_P10_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P20_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P25_2plus",
	"FIRST_RETURNS_elev_P30_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P40_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P50_2plus",
	"FIRST_RETURNS_elev_P60_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P70_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P75_2plus",
	"FIRST_RETURNS_elev_P80_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P90_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_elev_P95_2plus",
	"FIRST_RETURNS_elev_P99_2plus",
	"FIRST_RETURNS_elev_skewness_2plus",
	"FIRST_RETURNS_elev_stddev_2plus",
	"FIRST_RETURNS_elev_variance_2plus",
	"FIRST_RETURNS_int_AAD_2plus",
	"FIRST_RETURNS_int_ave_2plus",
	"FIRST_RETURNS_int_CV_2plus",
	"FIRST_RETURNS_int_IQ_2plus",
	"FIRST_RETURNS_int_kurtosis_2plus",
	"FIRST_RETURNS_int_max_2plus",
	"FIRST_RETURNS_int_min_2plus",
	"FIRST_RETURNS_int_mode_2plus",
	"FIRST_RETURNS_int_P01_2plus",
	"FIRST_RETURNS_int_P05_2plus",
	"FIRST_RETURNS_int_P10_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_int_P20_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_int_P25_2plus",
	"FIRST_RETURNS_int_P30_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_int_P40_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_int_P50_2plus",
	"FIRST_RETURNS_int_P60_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_int_P70_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_int_P75_2plus",
	"FIRST_RETURNS_int_P80_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_int_P90_2plus", # ADDED JUNE 2022
	"FIRST_RETURNS_int_P95_2plus",
	"FIRST_RETURNS_int_P99_2plus",
	"FIRST_RETURNS_int_skewness_2plus",
	"FIRST_RETURNS_int_stddev_2plus",
	"FIRST_RETURNS_int_variance_2plus",
	"Pct1stRtns_above_2",
	"Pct1stRtns_above_mean",
	"Pct1stRtns_above_mode", # ADDED JUNE 2022
	"PctAllRtns_above_2",
	"PctAllRtns_above_mean",
	"PctAllRtns_above_mode", # ADDED JUNE 2022
	"AllRtns_above_2_dividedby_Tot1stRtns_times100", # Added by TCB July 2022
    "AllRtns_above_mean_dividedby_Tot1stRtns_times100", # Added by TCB July 2022
    "AllRtns_above_mode_dividedby_Tot1stRtns_times100" # Added by TCB July 2022
	)

# ------- Strata metrics
StrataMetrics = function(z, r) {

	gt0 = which(z > 0)
	s0 = which(z > 0 & z < 0.5)
	s1 = which(z > 0.5 & z < 1)
	s2 = which(z > 1 & z < 2)
	s3 = which(z > 2 & z < 4)
	s4 = which(z > 4 & z < 8)
	s5 = which(z > 8 & z < 16)
	s6 = which(z > 16 & z < 32)
	s7 = which(z > 32 & z < 48)
	s8 = which(z > 48 & z < 64)
	s9 = which(z > 64)

	gt0first = which(z > 0 & r == 1)
	s0first = which(z > 0 & z < 0.5 & r == 1)
	s1first = which(z > 0.5 & z < 1 & r == 1)
	s2first = which(z > 1 & z < 2 & r == 1)
	s3first = which(z > 2 & z < 4 & r == 1)
	s4first = which(z > 4 & z < 8 & r == 1)
	s5first = which(z > 8 & z < 16 & r == 1)
	s6first = which(z > 16 & z < 32 & r == 1)
	s7first = which(z > 32 & z < 48 & r == 1)
	s8first = which(z > 48 & z < 64 & r == 1)
	s9first = which(z > 64 & r == 1)

  metrics = list(
	ALL_RETURNS_strata_0p5to1M_CV = sd(z[s1])/mean(z[s1]),
	ALL_RETURNS_strata_0p5to1M_kurtosis = kurtosis(z[s1]),
	ALL_RETURNS_strata_0p5to1M_return_proportion = length(z[s1]) / length(z[gt0]),
	ALL_RETURNS_strata_0p5to1M_skewness = skewness(z[s1]),
	ALL_RETURNS_strata_0p5to1M_stddev = sd(z[s1]),
	ALL_RETURNS_strata_0p5to1M_total_return_cnt = length(z[s1]),
	ALL_RETURNS_strata_0to0p5M_CV = sd(z[s0])/mean(z[s0]),
	ALL_RETURNS_strata_0to0p5M_kurtosis = kurtosis(z[s0]),
	ALL_RETURNS_strata_0to0p5M_return_proportion = length(z[s0]) / length(z[gt0]),
	ALL_RETURNS_strata_0to0p5M_skewness = skewness(z[s0]),
	ALL_RETURNS_strata_0to0p5M_stddev = sd(z[s0]),
	ALL_RETURNS_strata_0to0p5M_total_return_cnt = length(z[s0]),
	ALL_RETURNS_strata_16to32M_CV = sd(z[s6])/mean(z[s6]),
	ALL_RETURNS_strata_16to32M_kurtosis = kurtosis(z[s6]),
	ALL_RETURNS_strata_16to32M_return_proportion = length(z[s6]) / length(z[gt0]),
	ALL_RETURNS_strata_16to32M_skewness = skewness(z[s6]),
	ALL_RETURNS_strata_16to32M_stddev = sd(z[s6]),
	ALL_RETURNS_strata_16to32M_total_return_cnt = length(z[s6]),
	ALL_RETURNS_strata_1to2M_CV = sd(z[s2])/mean(z[s2]),
	ALL_RETURNS_strata_1to2M_kurtosis = kurtosis(z[s2]),
	ALL_RETURNS_strata_1to2M_return_proportion = length(z[s2]) / length(z[gt0]),
	ALL_RETURNS_strata_1to2M_skewness = skewness(z[s2]),
	ALL_RETURNS_strata_1to2M_stddev = sd(z[s2]),
	ALL_RETURNS_strata_1to2M_total_return_cnt = length(z[s2]),
	ALL_RETURNS_strata_2to4M_CV = sd(z[s3])/mean(z[s3]),
	ALL_RETURNS_strata_2to4M_kurtosis = kurtosis(z[s3]),
	ALL_RETURNS_strata_2to4M_return_proportion = length(z[s3]) / length(z[gt0]),
	ALL_RETURNS_strata_2to4M_skewness = skewness(z[s3]),
	ALL_RETURNS_strata_2to4M_stddev = sd(z[s3]),
	ALL_RETURNS_strata_2to4M_total_return_cnt = length(z[s3]),
	ALL_RETURNS_strata_32to48M_CV = sd(z[s7])/mean(z[s7]),
	ALL_RETURNS_strata_32to48M_kurtosis = kurtosis(z[s7]),
	ALL_RETURNS_strata_32to48M_return_proportion = length(z[s7]) / length(z[gt0]),
	ALL_RETURNS_strata_32to48M_skewness = skewness(z[s7]),
	ALL_RETURNS_strata_32to48M_stddev = sd(z[s7]),
	ALL_RETURNS_strata_32to48M_total_return_cnt = length(z[s7]),
	ALL_RETURNS_strata_48to64M_CV = sd(z[s8])/mean(z[s8]),
	ALL_RETURNS_strata_48to64M_kurtosis = kurtosis(z[s8]),
	ALL_RETURNS_strata_48to64M_return_proportion = length(z[s8]) / length(z[gt0]),
	ALL_RETURNS_strata_48to64M_skewness = skewness(z[s8]),
	ALL_RETURNS_strata_48to64M_stddev = sd(z[s8]),
	ALL_RETURNS_strata_48to64M_total_return_cnt = length(z[s8]),
	ALL_RETURNS_strata_4to8M_CV = sd(z[s4])/mean(z[s4]),
	ALL_RETURNS_strata_4to8M_kurtosis = kurtosis(z[s4]),
	ALL_RETURNS_strata_4to8M_return_proportion = length(z[s4]) / length(z[gt0]),
	ALL_RETURNS_strata_4to8M_skewness = skewness(z[s4]),
	ALL_RETURNS_strata_4to8M_stddev = sd(z[s4]),
	ALL_RETURNS_strata_4to8M_total_return_cnt = length(z[s4]),
	ALL_RETURNS_strata_64M_plus_CV = sd(z[s9])/mean(z[s9]),
	ALL_RETURNS_strata_64M_plus_kurtosis = kurtosis(z[s9]),
	ALL_RETURNS_strata_64M_plus_return_proportion = length(z[s9]) / length(z[gt0]),
	ALL_RETURNS_strata_64M_plus_skewness = skewness(z[s9]),
	ALL_RETURNS_strata_64M_plus_stddev = sd(z[s9]),
	ALL_RETURNS_strata_64M_plus_total_return_cnt = length(z[s9]),
	ALL_RETURNS_strata_8to16M_CV = sd(z[s5])/mean(z[s5]),
	ALL_RETURNS_strata_8to16M_kurtosis = kurtosis(z[s5]),
	ALL_RETURNS_strata_8to16M_return_proportion = length(z[s5]) / length(z[gt0]),
	ALL_RETURNS_strata_8to16M_skewness = skewness(z[s5]),
	ALL_RETURNS_strata_8to16M_stddev = sd(z[s5]),
	ALL_RETURNS_strata_8to16M_total_return_cnt = length(z[s5]),
	FIRST_RETURNS_strata_0p5to1M_CV = sd(z[s1first])/mean(z[s1first]),
	FIRST_RETURNS_strata_0p5to1M_kurtosis = kurtosis(z[s1first]),
	FIRST_RETURNS_strata_0p5to1M_return_proportion = length(z[s1first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_0p5to1M_skewness = skewness(z[s1first]),
	FIRST_RETURNS_strata_0p5to1M_stddev = sd(z[s1first]),
	FIRST_RETURNS_strata_0p5to1M_total_return_cnt = length(z[s1first]),
	FIRST_RETURNS_strata_0to0p5M_CV = sd(z[s0first])/mean(z[s0first]),
	FIRST_RETURNS_strata_0to0p5M_kurtosis = kurtosis(z[s0first]),
	FIRST_RETURNS_strata_0to0p5M_return_proportion = length(z[s0first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_0to0p5M_skewness = skewness(z[s0first]),
	FIRST_RETURNS_strata_0to0p5M_stddev = sd(z[s0first]),
	FIRST_RETURNS_strata_0to0p5M_total_return_cnt = length(z[s0first]),
	FIRST_RETURNS_strata_16to32M_CV = sd(z[s6first])/mean(z[s6first]),
	FIRST_RETURNS_strata_16to32M_kurtosis = kurtosis(z[s6first]),
	FIRST_RETURNS_strata_16to32M_return_proportion = length(z[s6first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_16to32M_skewness = skewness(z[s6first]),
	FIRST_RETURNS_strata_16to32M_stddev = sd(z[s6first]),
	FIRST_RETURNS_strata_16to32M_total_return_cnt = length(z[s6first]),
	FIRST_RETURNS_strata_1to2M_CV = sd(z[s2first])/mean(z[s2first]),
	FIRST_RETURNS_strata_1to2M_kurtosis = kurtosis(z[s2first]),
	FIRST_RETURNS_strata_1to2M_return_proportion = length(z[s2first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_1to2M_skewness = skewness(z[s2first]),
	FIRST_RETURNS_strata_1to2M_stddev = sd(z[s2first]),
	FIRST_RETURNS_strata_1to2M_total_return_cnt = length(z[s2first]),
	FIRST_RETURNS_strata_2to4M_CV = sd(z[s3first])/mean(z[s3first]),
	FIRST_RETURNS_strata_2to4M_kurtosis = kurtosis(z[s3first]),
	FIRST_RETURNS_strata_2to4M_return_proportion = length(z[s3first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_2to4M_skewness = skewness(z[s3first]),
	FIRST_RETURNS_strata_2to4M_stddev = sd(z[s3first]),
	FIRST_RETURNS_strata_2to4M_total_return_cnt = length(z[s3first]),
	FIRST_RETURNS_strata_32to48M_CV = sd(z[s7first])/mean(z[s7first]),
	FIRST_RETURNS_strata_32to48M_kurtosis = kurtosis(z[s7first]),
	FIRST_RETURNS_strata_32to48M_return_proportion = length(z[s7first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_32to48M_skewness = skewness(z[s7first]),
	FIRST_RETURNS_strata_32to48M_stddev = sd(z[s7first]),
	FIRST_RETURNS_strata_32to48M_total_return_cnt = length(z[s7first]),
	FIRST_RETURNS_strata_48to64M_CV = sd(z[s8first])/mean(z[s8first]),
	FIRST_RETURNS_strata_48to64M_kurtosis = kurtosis(z[s8first]),
	FIRST_RETURNS_strata_48to64M_return_proportion = length(z[s8first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_48to64M_skewness = skewness(z[s8first]),
	FIRST_RETURNS_strata_48to64M_stddev = sd(z[s8first]),
	FIRST_RETURNS_strata_48to64M_total_return_cnt = length(z[s8first]),
	FIRST_RETURNS_strata_4to8M_CV = sd(z[s4first])/mean(z[s4first]),
	FIRST_RETURNS_strata_4to8M_kurtosis = kurtosis(z[s4first]),
	FIRST_RETURNS_strata_4to8M_return_proportion = length(z[s4first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_4to8M_skewness = skewness(z[s4first]),
	FIRST_RETURNS_strata_4to8M_stddev = sd(z[s4first]),
	FIRST_RETURNS_strata_4to8M_total_return_cnt = length(z[s4first]),
	FIRST_RETURNS_strata_64M_plus_CV = sd(z[s9first])/mean(z[s9first]),
	FIRST_RETURNS_strata_64M_plus_kurtosis = kurtosis(z[s9first]),
	FIRST_RETURNS_strata_64M_plus_return_proportion = length(z[s9first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_64M_plus_skewness = skewness(z[s9first]),
	FIRST_RETURNS_strata_64M_plus_stddev = sd(z[s9first]),
	FIRST_RETURNS_strata_64M_plus_total_return_cnt = length(z[s9first]),
	FIRST_RETURNS_strata_8to16M_CV = sd(z[s5first])/mean(z[s5first]),
	FIRST_RETURNS_strata_8to16M_kurtosis = kurtosis(z[s5first]),
	FIRST_RETURNS_strata_8to16M_return_proportion = length(z[s5first]) / length(z[gt0first]),
	FIRST_RETURNS_strata_8to16M_skewness = skewness(z[s5first]),
	FIRST_RETURNS_strata_8to16M_stddev = sd(z[s5first]),
	FIRST_RETURNS_strata_8to16M_total_return_cnt = length(z[s5first])
	)

   return(metrics)
}

StrataMetricNames = c(
	"ALL_RETURNS_strata_0p5to1M_CV",
	"ALL_RETURNS_strata_0p5to1M_kurtosis",
	"ALL_RETURNS_strata_0p5to1M_return_proportion",
	"ALL_RETURNS_strata_0p5to1M_skewness",
	"ALL_RETURNS_strata_0p5to1M_stddev",
	"ALL_RETURNS_strata_0p5to1M_total_return_cnt",
	"ALL_RETURNS_strata_0to0p5M_CV",
	"ALL_RETURNS_strata_0to0p5M_kurtosis",
	"ALL_RETURNS_strata_0to0p5M_return_proportion",
	"ALL_RETURNS_strata_0to0p5M_skewness",
	"ALL_RETURNS_strata_0to0p5M_stddev",
	"ALL_RETURNS_strata_0to0p5M_total_return_cnt",
	"ALL_RETURNS_strata_16to32M_CV",
	"ALL_RETURNS_strata_16to32M_kurtosis",
	"ALL_RETURNS_strata_16to32M_return_proportion",
	"ALL_RETURNS_strata_16to32M_skewness",
	"ALL_RETURNS_strata_16to32M_stddev",
	"ALL_RETURNS_strata_16to32M_total_return_cnt",
	"ALL_RETURNS_strata_1to2M_CV",
	"ALL_RETURNS_strata_1to2M_kurtosis",
	"ALL_RETURNS_strata_1to2M_return_proportion",
	"ALL_RETURNS_strata_1to2M_skewness",
	"ALL_RETURNS_strata_1to2M_stddev",
	"ALL_RETURNS_strata_1to2M_total_return_cnt",
	"ALL_RETURNS_strata_2to4M_CV",
	"ALL_RETURNS_strata_2to4M_kurtosis",
	"ALL_RETURNS_strata_2to4M_return_proportion",
	"ALL_RETURNS_strata_2to4M_skewness",
	"ALL_RETURNS_strata_2to4M_stddev",
	"ALL_RETURNS_strata_2to4M_total_return_cnt",
	"ALL_RETURNS_strata_32to48M_CV",
	"ALL_RETURNS_strata_32to48M_kurtosis",
	"ALL_RETURNS_strata_32to48M_return_proportion",
	"ALL_RETURNS_strata_32to48M_skewness",
	"ALL_RETURNS_strata_32to48M_stddev",
	"ALL_RETURNS_strata_32to48M_total_return_cnt",
	"ALL_RETURNS_strata_48to64M_CV",
	"ALL_RETURNS_strata_48to64M_kurtosis",
	"ALL_RETURNS_strata_48to64M_return_proportion",
	"ALL_RETURNS_strata_48to64M_skewness",
	"ALL_RETURNS_strata_48to64M_stddev",
	"ALL_RETURNS_strata_48to64M_total_return_cnt",
	"ALL_RETURNS_strata_4to8M_CV",
	"ALL_RETURNS_strata_4to8M_kurtosis",
	"ALL_RETURNS_strata_4to8M_return_proportion",
	"ALL_RETURNS_strata_4to8M_skewness",
	"ALL_RETURNS_strata_4to8M_stddev",
	"ALL_RETURNS_strata_4to8M_total_return_cnt",
	"ALL_RETURNS_strata_64M_plus_CV",
	"ALL_RETURNS_strata_64M_plus_kurtosis",
	"ALL_RETURNS_strata_64M_plus_return_proportion",
	"ALL_RETURNS_strata_64M_plus_skewness",
	"ALL_RETURNS_strata_64M_plus_stddev",
	"ALL_RETURNS_strata_64M_plus_total_return_cnt",
	"ALL_RETURNS_strata_8to16M_CV",
	"ALL_RETURNS_strata_8to16M_kurtosis",
	"ALL_RETURNS_strata_8to16M_return_proportion",
	"ALL_RETURNS_strata_8to16M_skewness",
	"ALL_RETURNS_strata_8to16M_stddev",
	"ALL_RETURNS_strata_8to16M_total_return_cnt",
	"FIRST_RETURNS_strata_0p5to1M_CV",
	"FIRST_RETURNS_strata_0p5to1M_kurtosis",
	"FIRST_RETURNS_strata_0p5to1M_return_proportion",
	"FIRST_RETURNS_strata_0p5to1M_skewness",
	"FIRST_RETURNS_strata_0p5to1M_stddev",
	"FIRST_RETURNS_strata_0p5to1M_total_return_cnt",
	"FIRST_RETURNS_strata_0to0p5M_CV",
	"FIRST_RETURNS_strata_0to0p5M_kurtosis",
	"FIRST_RETURNS_strata_0to0p5M_return_proportion",
	"FIRST_RETURNS_strata_0to0p5M_skewness",
	"FIRST_RETURNS_strata_0to0p5M_stddev",
	"FIRST_RETURNS_strata_0to0p5M_total_return_cnt",
	"FIRST_RETURNS_strata_16to32M_CV",
	"FIRST_RETURNS_strata_16to32M_kurtosis",
	"FIRST_RETURNS_strata_16to32M_return_proportion",
	"FIRST_RETURNS_strata_16to32M_skewness",
	"FIRST_RETURNS_strata_16to32M_stddev",
	"FIRST_RETURNS_strata_16to32M_total_return_cnt",
	"FIRST_RETURNS_strata_1to2M_CV",
	"FIRST_RETURNS_strata_1to2M_kurtosis",
	"FIRST_RETURNS_strata_1to2M_return_proportion",
	"FIRST_RETURNS_strata_1to2M_skewness",
	"FIRST_RETURNS_strata_1to2M_stddev",
	"FIRST_RETURNS_strata_1to2M_total_return_cnt",
	"FIRST_RETURNS_strata_2to4M_CV",
	"FIRST_RETURNS_strata_2to4M_kurtosis",
	"FIRST_RETURNS_strata_2to4M_return_proportion",
	"FIRST_RETURNS_strata_2to4M_skewness",
	"FIRST_RETURNS_strata_2to4M_stddev",
	"FIRST_RETURNS_strata_2to4M_total_return_cnt",
	"FIRST_RETURNS_strata_32to48M_CV",
	"FIRST_RETURNS_strata_32to48M_kurtosis",
	"FIRST_RETURNS_strata_32to48M_return_proportion",
	"FIRST_RETURNS_strata_32to48M_skewness",
	"FIRST_RETURNS_strata_32to48M_stddev",
	"FIRST_RETURNS_strata_32to48M_total_return_cnt",
	"FIRST_RETURNS_strata_48to64M_CV",
	"FIRST_RETURNS_strata_48to64M_kurtosis",
	"FIRST_RETURNS_strata_48to64M_return_proportion",
	"FIRST_RETURNS_strata_48to64M_skewness",
	"FIRST_RETURNS_strata_48to64M_stddev",
	"FIRST_RETURNS_strata_48to64M_total_return_cnt",
	"FIRST_RETURNS_strata_4to8M_CV",
	"FIRST_RETURNS_strata_4to8M_kurtosis",
	"FIRST_RETURNS_strata_4to8M_return_proportion",
	"FIRST_RETURNS_strata_4to8M_skewness",
	"FIRST_RETURNS_strata_4to8M_stddev",
	"FIRST_RETURNS_strata_4to8M_total_return_cnt",
	"FIRST_RETURNS_strata_64M_plus_CV",
	"FIRST_RETURNS_strata_64M_plus_kurtosis",
	"FIRST_RETURNS_strata_64M_plus_return_proportion",
	"FIRST_RETURNS_strata_64M_plus_skewness",
	"FIRST_RETURNS_strata_64M_plus_stddev",
	"FIRST_RETURNS_strata_64M_plus_total_return_cnt",
	"FIRST_RETURNS_strata_8to16M_CV",
	"FIRST_RETURNS_strata_8to16M_kurtosis",
	"FIRST_RETURNS_strata_8to16M_return_proportion",
	"FIRST_RETURNS_strata_8to16M_skewness",
	"FIRST_RETURNS_strata_8to16M_stddev",
	"FIRST_RETURNS_strata_8to16M_total_return_cnt"
	)

# ------- Canopy height model-based metrics
# Adapted from https://jean-romain.github.io/lidRbook/outbox.html#outbox-custom-metrics

# Rumple
grid_rumple_index <- function(las, res, start) {
  las <- filter_surfacepoints(las, 1)
  return(grid_metrics(las, ~rumple_index(X, Y, Z), res, start))
}

# Mean CHM
grid_chm_mean <- function(las, res, start) {
  las <- filter_surfacepoints(las, 1)
  return(grid_metrics(las, ~mean(Z), res, start))
}

# SD CHM
grid_chm_sd <- function(las, res, start) {
  las <- filter_surfacepoints(las, 1)
  return(grid_metrics(las, ~sd(Z), res, start))
}

# Max CHM
grid_chm_max <- function(las, res, start) {
  las <- filter_surfacepoints(las, 1)
  return(grid_metrics(las, ~max(Z), res, start))
}

# CHM FVP
# FUSION's surface_volume_ratio of GridSurfaceStats (p. 94 of FUSION manual)
# surface volume ratio = surface volume / potential volume
# here I have defined surface volume as sum(Z), where Z are 1-m CHM heights (note that the CHM must be 1-m), and
# I have defined potential volume as res*res*max(Z)
# FUSION uses a more refined TIN surface for surface volume, but I don't know how to code this within lidR
# NOTE: THIS DOESN'T WORK AS INTENDED!!!!
grid_chm_fpv <- function(las, res, start) {
  las <- filter_surfacepoints(las, 1)
  return(grid_metrics(las, ~sum(Z)/(res*res*max(Z)), res, start))
}

# Canopy height model metrics routine
CHMmetrics <- function(chunk, res, start, M) {
  bbox <- raster::extent(chunk)
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  if (M == "rumple") metric <- grid_rumple_index(las, res, start)
  if (M == "chm_mean") metric <- grid_chm_mean(las, res, start)
  if (M == "chm_sd") metric <- grid_chm_sd(las, res, start)
  if (M == "chm_max") metric <- grid_chm_max(las, res, start)
  if (M == "chm_fpv") metric <- grid_chm_fpv(las, res, start)
  metric <- raster::crop(metric, bbox)
  return(metric)
}

CHMMetricNames = c("canopy_average_height","canopy_FPV","canopy_maximum_height","canopy_rumple","canopy_stddev_height")
