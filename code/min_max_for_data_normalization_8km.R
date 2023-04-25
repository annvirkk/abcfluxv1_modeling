



### Prepping

.libPaths("/home/master/R/x86_64-pc-linux-gnu-library/4.2")


library("dplyr")
library("caret")
library("parallel")
library("doParallel")
library("randomForest")
library("quantregForest")  
library("groupdata2")  
library(itertools)
library(foreach)
library(doParallel)
library(parallel)
library(kernlab)
library(nnet)
library("sp")
library("ggplot2")
library("caret")
library("dplyr")
library("purrr")
library("raster")
library("terra")
library(stringr)
library(googleCloudStorageR)
library(purrr)


# Terra settings
terraOptions(memfrac=0.9, tempdir = "/home/master/temp/") 
options(digits=20) # this is important, otherwise we will lose some digits with the cropped extents


# Google cloud settings
my_project_id <- "top-operand-328213"
gcs_list_buckets(my_project_id)
gcs_global_bucket("abcflux_modeling_files")
contents <- gcs_list_objects()
gcs_upload_set_limit(50000000L) # increasing data size limit for transferring data to google cloud




### Load the model training data
d <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg_allsites.csv")
d <- d %>% mutate_if(is.integer, as.numeric)



# Empty data frame for storing
df <- data.frame(nrow=1)



### Normalize model training data so that extreme values would not get as much emphasis.
# A common procedure in ML works
# Nomalization between 0 and 1
# We need to consider the range of values both at the sites (model training data) as well as across the entire domain (predictors)


# Load all the monthly predictors from 1882-2016 for this purpose and get their max and min

gcs_get_object("predictors_8km/soc.tif", saveToDisk = "predictors_8km/soc.tif", overwrite=TRUE)
SoilGrids_SOC_SoilGrids_SOCstock <-  rast("predictors_8km/soc.tif")
#plot(SoilGrids_SOC_SoilGrids_SOCstock)
SoilGrids_SOC_SoilGrids_SOCstock
summary(d$SoilGrids_SOC_SoilGrids_SOCstock)
SoilGrids_SOC_SoilGrids_SOCstock <- SoilGrids_SOC_SoilGrids_SOCstock/100
#SoilGrids_SOC_SoilGrids_SOCstock <- as.data.frame(SoilGrids_SOC_SoilGrids_SOCstock, xy=TRUE)
# Unit tonnes per ha
socmin <- min(c(d$SoilGrids_SOC_SoilGrids_SOCstock, values(SoilGrids_SOC_SoilGrids_SOCstock)), na.rm=TRUE)
socmax <- max(c(d$SoilGrids_SOC_SoilGrids_SOCstock, values(SoilGrids_SOC_SoilGrids_SOCstock)), na.rm=TRUE)
df <- cbind(df, socmin)
df <- cbind(df, socmax)

gcs_get_object("predictors_8km/cti.tif", saveToDisk = "predictors_8km/cti.tif", overwrite=TRUE)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- rast("predictors_8km/cti.tif")
#plot(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m
summary(d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m/100

#dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- as.data.frame(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m, xy=TRUE)
# Compound topographic index, high value= high topographic moisture
dtmmin <- min(c(d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m, values(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)), na.rm=TRUE)
dtmmax <- max(c(d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m, values(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)), na.rm=TRUE)
df <- cbind(df, dtmmin)
df <- cbind(df, dtmmax)

gcs_get_object("predictors_8km/ph.tif", saveToDisk = "predictors_8km/ph.tif", overwrite=TRUE)
PHIHOX_M_sl1_250m_ll_SoilGrids <- rast("predictors_8km/ph.tif")
#plot(PHIHOX_M_sl1_250m_ll_SoilGrids)
PHIHOX_M_sl1_250m_ll_SoilGrids
summary(d$PHIHOX_M_sl1_250m_ll_SoilGrids)
PHIHOX_M_sl1_250m_ll_SoilGrids <- PHIHOX_M_sl1_250m_ll_SoilGrids/100
#PHIHOX_M_sl1_250m_ll_SoilGrids <- as.data.frame(PHIHOX_M_sl1_250m_ll_SoilGrids, xy=TRUE)
phmin <- min(c(d$PHIHOX_M_sl1_250m_ll_SoilGrids, values(PHIHOX_M_sl1_250m_ll_SoilGrids)), na.rm=TRUE)
phmax <- max(c(d$PHIHOX_M_sl1_250m_ll_SoilGrids, values(PHIHOX_M_sl1_250m_ll_SoilGrids)), na.rm=TRUE)
df <- cbind(df, phmin)
df <- cbind(df, phmax)

gcs_get_object("predictors_8km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif", saveToDisk = "predictors_8km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif", overwrite=TRUE)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- rast("predictors_8km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif")
#plot(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH
summary(d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH/100
#UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- as.data.frame(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH, xy=TRUE)
# Permafrost probability (fraction 0-1)
permamin <- min(c(d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH, values(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)), na.rm=TRUE)
permamax <- max(c(d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH, values(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)), na.rm=TRUE)
df <- cbind(df, permamin)
df <- cbind(df, permamax)


# Dynamic rasters
# Load all the data from the different years, 1982 to 2016
# Need to check that the model training data and prediction data are the same

# List all the data in Google bucket
contents <- gcs_list_objects()
files_to_download <- grep("*.tif", contents$name, value = TRUE)

# working directory 
setwd("/home/master/")


# Create a function for downloading the data
select_data <- function(var) {
  
  files_to_download2 <- files_to_download[grepl(var, files_to_download)]
  files_to_download3 <- files_to_download2[!(grepl(paste(var, "195", sep="_"), files_to_download2)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "196", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "197", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2017", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2018", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2019", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2020", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2021", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "1980", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl("trend", files_to_download3)) ]
  files_to_download3 <- files_to_download3[endsWith(files_to_download3, '.tif')] 
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "1981", sep="_"), files_to_download3)) ] 
  files_to_download3 <<- files_to_download3 # <<- saves the output
  
  map(files_to_download3, function(x) gcs_get_object(x, saveToDisk = x, overwrite = TRUE))
  
  
  
}



# srad
select_data(var="predictors_8km/srad")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE); minr[,1]/10
maxr <- global(r, "max", na.rm=TRUE) ; maxr[,1]/10
# hist(values(rast(files_to_download3[1])/10))
# test <- d %>% filter(Interval==1); hist(test$srad_terraclimate_sites)
sradmin <- min(c(d$srad_terraclimate_sites, minr[,1]), na.rm=TRUE)
sradmax <- max(c(d$srad_terraclimate_sites, maxr[,1]), na.rm=TRUE)
df <- cbind(df, sradmin)
df <- cbind(df, sradmax)




# tmean
select_data(var="predictors_8km/tmean")
print(files_to_download3)
r <- rast(files_to_download3)
# hist(values(rast(files_to_download3[1])/10))
# test <- d %>% filter(Interval==1); hist(test$tmean_terraclimate_sites)
minr <- global(r, "min", na.rm=TRUE); minr[,1]/10
maxr <- global(r, "max", na.rm=TRUE) ; maxr[,1]/10
tmeanmin <- min(c(d$tmean_terraclimate_sites, minr[,1]), na.rm=TRUE)
tmeanmax <- max(c(d$tmean_terraclimate_sites, maxr[,1]), na.rm=TRUE)
df <- cbind(df, tmeanmin)
df <- cbind(df, tmeanmax)

# vpd
select_data(var="predictors_8km/vpd")
print(files_to_download3)
r <- rast(files_to_download3)
# hist(values(rast(files_to_download3[1])))
# test <- d %>% filter(Interval==1); hist(test$vpd_terraclimate_sites)
minr <- global(r, "min", na.rm=TRUE)
maxr <- global(r, "max", na.rm=TRUE) 
vpdmin <- min(c(d$vpd_terraclimate_sites, minr[,1]), na.rm=TRUE)
vpdmax <- max(c(d$vpd_terraclimate_sites, maxr[,1]), na.rm=TRUE)
df <- cbind(df, vpdmin)
df <- cbind(df, vpdmax)



select_data(var="predictors_8km/co2")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE); minr[,1]/1000
maxr <- global(r, "max", na.rm=TRUE) ; maxr[,1]/1000
# hist(values(rast(files_to_download3[1]))/1000 )
# test <- d %>% filter(Interval==1); hist(test$Barrow_CO2_conc_Barrow_CO2conc)
co2min <- min(c(d$Barrow_CO2_conc_Barrow_CO2conc , minr[,1]), na.rm=TRUE)
co2max <- max(c(d$Barrow_CO2_conc_Barrow_CO2conc , maxr[,1]), na.rm=TRUE)
df <- cbind(df, co2min)
df <- cbind(df, co2max)

select_data(var="predictors_8km/ndvi_gimms")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE); minr[,1]/10000
maxr <- global(r, "max", na.rm=TRUE) ; maxr[,1]/10000
# hist(values(rast(files_to_download3[1]))/10000 )
# test <- d %>% filter(Interval==1); hist(test$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled)
ndvimin <- min(c(d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled , minr[,1]), na.rm=TRUE)
ndvimax <- max(c(d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled , maxr[,1]), na.rm=TRUE)
df <- cbind(df, ndvimin)
df <- cbind(df, ndvimax)

select_data(var="predictors_8km/soiltemplevel1")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE); minr[,1]/100
maxr <- global(r, "max", na.rm=TRUE) ; maxr[,1]/100
# hist(values(rast(files_to_download3[1]))/100 )
# test <- d %>% filter(Interval==1); hist(test$Soil.temperature.level.1_era5_soilmoist_temp_snow)
soiltempmin <- min(c(d$Soil.temperature.level.1_era5_soilmoist_temp_snow , minr[,1]), na.rm=TRUE)
soiltempmax <- max(c(d$Soil.temperature.level.1_era5_soilmoist_temp_snow , maxr[,1]), na.rm=TRUE)
df <- cbind(df, soiltempmin)
df <- cbind(df, soiltempmax)



select_data(var="predictors_8km/snowcover")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE); minr[,1]/100
maxr <- global(r, "max", na.rm=TRUE) ; maxr[,1]/100
# hist(values(rast(files_to_download3[1]))/100 )
# test <- d %>% filter(Interval==1); hist(test$Snow.cover_era5_soilmoist_temp_snow)
snowcmin <- min(c(d$Snow.cover_era5_soilmoist_temp_snow , minr[, 1]), na.rm=TRUE)
snowcmax <- max(c(d$Snow.cover_era5_soilmoist_temp_snow , maxr[, 1]), na.rm=TRUE)
df <- cbind(df, snowcmin)
df <- cbind(df, snowcmax)


select_data(var="predictors_8km/soilmoistlevel1")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE); minr[,1]/100
maxr <- global(r, "max", na.rm=TRUE) ; maxr[,1]/100
# hist(values(rast(files_to_download3[1]))/100 )
# test <- d %>% filter(Interval==1); hist(test$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
soilmoistmin <- min(c(d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow , minr[,1]), na.rm=TRUE)
soilmoistmax <- max(c(d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow , maxr[,1]), na.rm=TRUE)
df <- cbind(df, soilmoistmin)
df <- cbind(df, soilmoistmax)




select_data(var="predictors_8km/Percent_TreeCover_VCF5KYR")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE)
maxr <- global(r, "max", na.rm=TRUE) 
# hist(values(rast(files_to_download3[1])))
# test <- d %>% filter(Interval==1); hist(test$Percent_TreeCover_AVHRR_VCF5KYR)
treecmin <- min(c(d$Percent_TreeCover_AVHRR_VCF5KYR , minr[,1]), na.rm=TRUE)
treecmax <- max(c(d$Percent_TreeCover_AVHRR_VCF5KYR , maxr[,1]), na.rm=TRUE)
df <- cbind(df, treecmin)
df <- cbind(df, treecmax)

select_data(var="predictors_8km/Percent_NonTree_Vegetation_VCF5KYR")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE)
maxr <- global(r, "max", na.rm=TRUE)
# hist(values(rast(files_to_download3[1])))
# test <- d %>% filter(Interval==1); hist(test$Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
nontreecmin <- min(c(d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR  , minr[,1]), na.rm=TRUE)
nontreecmax <- max(c(d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR , maxr[,1]), na.rm=TRUE)
df <- cbind(df, nontreecmin)
df <- cbind(df, nontreecmax)


select_data(var="predictors_8km/Percent_NonVegetated_VCF5KYR")
print(files_to_download3)
r <- rast(files_to_download3)
minr <- global(r, "min", na.rm=TRUE)
maxr <- global(r, "max", na.rm=TRUE)
# hist(values(rast(files_to_download3[1])))
# test <- d %>% filter(Interval==1); hist(test$Percent_NonVegetated_AVHRR_VCF5KYR)
nonvegcmin <- min(c(d$Percent_NonVegetated_AVHRR_VCF5KYR  , minr[,1]), na.rm=TRUE)
nonvegcmax <- max(c(d$Percent_NonVegetated_AVHRR_VCF5KYR , maxr[,1]), na.rm=TRUE)
df <- cbind(df, nonvegcmin)
df <- cbind(df, nonvegcmax)


# write out
write.csv(df, "/home/master/abcfluxv1_modeling/results/minmax_8km.csv", row.names=FALSE)

