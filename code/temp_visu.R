

library("terra")

setwd("/home/master/predictions_8km/raster/0.5")


gpp <- rast("GPP_gC_m2_8km_1990_2016_sum.tif")

reco <- rast("Reco_gC_m2_8km_1990_2016_sum.tif")

nee <- reco-gpp

writeRaster(nee, "NEE_gC_m2_8km_1990_2016_sum_fromGPPReco.tif", overwrite=TRUE)


nee_orig <- rast("NEE_gC_m2_8km_1990_2016_sum.tif")

rr <- nee - nee_orig


# trend?
gpp <- rast("GPP_gC_m2_8km_1990_2016_sum_trend.tif")

reco <- rast("Reco_gC_m2_8km_1990_2016_sum_trend.tif")

nee_trend <- reco[[2]] - gpp[[2]]


writeRaster(nee_trend, "NEE_gC_m2_8km_1990_2016_sum_trend_fromGPPReco.tif", overwrite=TRUE)



setwd("/home/master/trends_8km/")

r <- rast("tmean_1990_2016_annual_trend.tif")
plot(r[[2]])

r <- rast("ndvi_1990_2016_annual_trend.tif")
plot(r[[2]])

r <- rast("vpd_1990_2016_annual_trend.tif")
plot(r[[2]])


