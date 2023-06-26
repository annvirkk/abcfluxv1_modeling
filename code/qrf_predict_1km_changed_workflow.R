


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
d <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg_allsites_1km.csv")
d <- d %>% mutate_if(is.integer, as.numeric)

# Reverse gpp first so that GPP and Reco have both positive values
d$GPP_gC_m2 <- -d$GPP_gC_m2 

# Simplifying the land cover classes (for simplicity and because some classes have very limited amount of data)
d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==1, 31, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)

d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==90, 30, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)

d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==41, 160, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)



# # Define the factors
# d$Study_ID_Short <- factor(d$Study_ID_Short)
# d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- as.factor(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
# d$Interval <- as.factor(d$Interval)



# ### Normalize model training data so that extreme values would not get as much emphasis.
# 
# ### Do the actual normalization  (value – min) / (max – min) 
# normalize <- function(x, ...) {
#   return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
# }
# d$NEE_gC_m2 <- normalize(d$NEE_gC_m2, na.rm=TRUE)
# d$Reco_gC_m2 <- normalize(d$Reco_gC_m2, na.rm=TRUE)
# d$GPP_gC_m2 <- normalize(d$GPP_gC_m2, na.rm=TRUE)


### Model parameters

# Predictors
Baseline_vars_1km <- c("srad_terraclimate_sites", "vpd_terraclimate_sites",
                       
                       "Barrow_CO2_conc_Barrow_CO2conc", # atmos CO2 conc
                       
                       "Snow.cover_era5_soilmoist_temp_snow",  "Soil.temperature.level.1_era5_soilmoist_temp_snow", "Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow", #era5 here
                       
                       "NDVI_whittaker_constant_monthly_mean", # Optical RS, dropped several because highly correlated
                       
                       "LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality", # surface temp
                       
                       
                       "dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m", 
                       
                       
                       "Percent_NonTree_Vegetation_MOD44B_sites", "Percent_NonVegetated_MOD44B_sites", "Percent_Tree_Cover_MOD44B_sites", # veg cover
                       
                       
                       "PHIHOX_M_sl1_250m_ll_SoilGrids", "SoilGrids_SOC_SoilGrids_SOCstock",
                       
                       "UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH" # Permafrost (static)
                       
                       
                       
)

# check that the columns exist
Baseline_vars_1km %in% colnames(d)



### Response variables
resp_vars <- c("NEE_gC_m2",  "GPP_gC_m2", "Reco_gC_m2") 




### Time periods for the monthly predictions
# Loop through average monthly rasters
time <- seq(as.Date("2001/01/01"), as.Date("2020/12/31"), "months")
time <- substr(time, 1, 7)
time <- sub("-", "_", sub("_", "", time, fixed=TRUE), fixed=TRUE)
time_alt <- gsub("_0", "_", time)

# 
# ## TEMPORARY
# time <- time[7]
# time_alt <- time_alt[7]

### Set up clusters for parallel processing
#  https://stackoverflow.com/questions/44774516/parallel-processing-in-r-in-caret
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
library(doParallel)
# create the cluster for caret to use
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/master/R/x86_64-pc-linux-gnu-library/4.2"))




# Quantiles to get for predictions
quantiles = c(0.025, 0.5, 0.975)



### Load static rasters and normalize

### Load static vars (only once)
setwd("/home/master/")


gcs_get_object("masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_northpolelambert1km_tundraboreal_simplevers.tif", saveToDisk = "masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_northpolelambert1km_tundraboreal_simplevers.tif",overwrite=TRUE)
ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- rast("masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_northpolelambert1km_tundraboreal_simplevers.tif")
# # reclassifying was done in Arcgis
# unique(values(ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged))
# 
# m <- c(1, 1, 31,
#        90, 90, 30,
#        41, 41, 160)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# rcx2 <- classify(ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged, rclmat, include.lowest=TRUE, others=NULL)
# ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- rcx2
# rm(rcx2); gc()
# # plot(rcx2)
# # hist(values(rcx2))
# # unique(rcx2)
# # 
# # identical to this, so correct
# m <- c(1, 1, 31,
#        21, 21, 21,
#        30, 30, 30,
#        31, 31, 31,
#        33, 33, 33,
#        60, 60, 60,
#        70, 70, 70,
#        80, 80, 80,
#        120, 120, 120,
#        90, 90, 30,
#        41, 41, 160, 
#        160, 160, 160)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# rcx2 <- classify(ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged, rclmat, include.lowest=TRUE)
# hist(values(rcx2))
# unique(values(rcx2))
# 
# 
# d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==1, 31, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
# 
# d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==90, 30, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
# 
# d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==41, 160, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)


gcs_get_object("predictors_1km/soc.tif", saveToDisk = "predictors_1km/soc.tif", overwrite=TRUE)
SoilGrids_SOC_SoilGrids_SOCstock <-  rast("predictors_1km/soc.tif")
#plot(SoilGrids_SOC_SoilGrids_SOCstock)
SoilGrids_SOC_SoilGrids_SOCstock
summary(d$SoilGrids_SOC_SoilGrids_SOCstock)
NAflag(SoilGrids_SOC_SoilGrids_SOCstock) <- 0 # checked from the original layer that there are no zero pixels
SoilGrids_SOC_SoilGrids_SOCstock <- SoilGrids_SOC_SoilGrids_SOCstock/100
#SoilGrids_SOC_SoilGrids_SOCstock <- as.data.frame(SoilGrids_SOC_SoilGrids_SOCstock, xy=TRUE)
# Unit tonnes per ha

gcs_get_object("predictors_1km/cti.tif", saveToDisk = "predictors_1km/cti.tif", overwrite=TRUE)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- rast("predictors_1km/cti.tif")
#plot(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m
summary(d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m/100
#dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- as.data.frame(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m, xy=TRUE)
# Compound topographic index, high value= high topographic moisture



gcs_get_object("predictors_1km/ph.tif", saveToDisk = "predictors_1km/ph.tif", overwrite=TRUE)
PHIHOX_M_sl1_250m_ll_SoilGrids <- rast("predictors_1km/ph.tif")
#plot(PHIHOX_M_sl1_250m_ll_SoilGrids)
PHIHOX_M_sl1_250m_ll_SoilGrids
summary(d$PHIHOX_M_sl1_250m_ll_SoilGrids)
PHIHOX_M_sl1_250m_ll_SoilGrids[PHIHOX_M_sl1_250m_ll_SoilGrids<= 0] <- NA
PHIHOX_M_sl1_250m_ll_SoilGrids <- PHIHOX_M_sl1_250m_ll_SoilGrids/100
#PHIHOX_M_sl1_250m_ll_SoilGrids <- as.data.frame(PHIHOX_M_sl1_250m_ll_SoilGrids, xy=TRUE)



gcs_get_object("predictors_1km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif", saveToDisk = "predictors_1km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif", overwrite=TRUE)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- rast("predictors_1km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif")
#plot(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH
summary(d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- mask(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH, ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged) # this only worked
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH/100
#UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- as.data.frame(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH, xy=TRUE)
# Permafrost probability (fraction 0-1)


# raster merging
pred_rast_static <- c( ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged, SoilGrids_SOC_SoilGrids_SOCstock,
                       dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m,
                       PHIHOX_M_sl1_250m_ll_SoilGrids,
                       UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH) 

# change raster names so that they are exactly the same as in the model file
names(pred_rast_static) <- c( "ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged", "SoilGrids_SOC_SoilGrids_SOCstock",
                              "dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m",
                              "PHIHOX_M_sl1_250m_ll_SoilGrids",
                              "UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH") 

rm(SoilGrids_SOC_SoilGrids_SOCstock)
rm(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)
rm(PHIHOX_M_sl1_250m_ll_SoilGrids)
rm(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)
gc()


# # convert to a data frame - no, this takes too much time
# pred_rast_static_df <- as.data.frame(pred_rast_static, xy=TRUE, na.rm=TRUE) # na.rm=TRUE removes NA rows
# 
# # remove if any of the variables have NA in the pixel
# pred_rast_static_na <- na.omit(pred_rast_static_df)
# str(pred_rast_static_na)



# ## now also make an additional fake raster to merge the full pred raster with that has all the model training data
# fake_rast <- crop(pred_rast_static, ext(0.5e+06, 0.5e+06+(500*999.97108203905224855), -3.5e+06, -3.5e+06+(500*999.97108203905224855)))
# ext(fake_rast) <- c(0.5e+06, 0.5e+06+(500*999.97108203905224855), -4.5e+06, -4.5e+06+(500*999.97108203905224855))
# #fake_rast <- resample(fake_rast, pred_rast_static, method="near")
# fake_rast[[1]] <- d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged # change esa cci veg type values
# #fake_rast[[5]] <- d$Interval # change interval values
# test <- merge(pred_rast_static, fake_rast) # can be heavy to compute see https://github.com/rspatial/terra/issues/210
# rm(fake_rast)


print("static vars loaded")



# loop through the time periods and load dynamic data rasters
for (t in 1:length(time)) {
  
  
  
  
  # t <- 7  
  setwd("/home/master/") 
  gcs_get_object(paste0("predictors_1km/", "srad_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_1km/", "srad_", time_alt[t], ".tif"), overwrite=TRUE)
  srad_terraclimate_sites <- rast(paste0("/home/master/predictors_1km/", "srad_", time_alt[t], ".tif"))
  #plot(srad_terraclimate_sites)
  srad_terraclimate_sites
  summary(d$srad_terraclimate_sites)
  srad_terraclimate_sites <- srad_terraclimate_sites/10  
  names(srad_terraclimate_sites) <- "srad_terraclimate_sites"
 #srad_terraclimate_sites <- as.data.frame(srad_terraclimate_sites, xy=TRUE)
  # Downward surface shortwave radiation. Unit W/m2. Both need to be divided by 10 to get to the original scale
  

  pred_rasters <- c(pred_rast_static, srad_terraclimate_sites)
  rm(srad_terraclimate_sites)
  rm(pred_rast_static)
  gc()
  
  gcs_get_object(paste0("predictors_1km/", "co2_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_1km/", "co2_", time_alt[t], ".tif"), overwrite=TRUE)
  Barrow_CO2_conc_Barrow_CO2conc <- rast(paste0("/home/master/predictors_1km/", "co2_", time_alt[t], ".tif"))
  #plot(Barrow_CO2_conc_Barrow_CO2conc)
  Barrow_CO2_conc_Barrow_CO2conc
  summary(d$Barrow_CO2_conc_Barrow_CO2conc)
  Barrow_CO2_conc_Barrow_CO2conc <- Barrow_CO2_conc_Barrow_CO2conc/1000  
  names(Barrow_CO2_conc_Barrow_CO2conc) <- "Barrow_CO2_conc_Barrow_CO2conc"
 #Barrow_CO2_conc_Barrow_CO2conc <- as.data.frame(Barrow_CO2_conc_Barrow_CO2conc, xy=TRUE)
  # atm CO2 concentrations in ppm
  
  pred_rasters <- c(pred_rasters, Barrow_CO2_conc_Barrow_CO2conc)
  rm(Barrow_CO2_conc_Barrow_CO2conc)
  gc()
  
  
  gcs_get_object(paste0("predictors_1km/", "ndvi_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_1km/", "ndvi_", time_alt[t], ".tif"), overwrite=TRUE)
  NDVI_whittaker_constant_monthly_mean <- rast(paste0("/home/master/predictors_1km/","ndvi_", time_alt[t], ".tif")) 
  #plot(NDVI_whittaker_constant_monthly_mean)
  NDVI_whittaker_constant_monthly_mean[NDVI_whittaker_constant_monthly_mean< -10000] <- NA
  NDVI_whittaker_constant_monthly_mean
  summary(d$NDVI_whittaker_constant_monthly_mean)
  NDVI_whittaker_constant_monthly_mean <- NDVI_whittaker_constant_monthly_mean/10000 
  names(NDVI_whittaker_constant_monthly_mean) <- "NDVI_whittaker_constant_monthly_mean"
  #NDVI_whittaker_constant_monthly_mean <- as.data.frame(NDVI_whittaker_constant_monthly_mean, xy=TRUE)
  
  pred_rasters <- c(pred_rasters, NDVI_whittaker_constant_monthly_mean)
  rm(NDVI_whittaker_constant_monthly_mean)
  gc()
  
  gcs_get_object(paste0("predictors_1km/", "soiltemplevel1_", time[t], ".tif"), saveToDisk = paste0("predictors_1km/", "soiltemplevel1_", time[t], ".tif"), overwrite=TRUE)
  Soil.temperature.level.1_era5_soilmoist_temp_snow <- rast(paste0("/home/master/predictors_1km/","soiltemplevel1_", time[t], ".tif"))
  #plot(Soil.temperature.level.1_era5_soilmoist_temp_snow)
  Soil.temperature.level.1_era5_soilmoist_temp_snow[Soil.temperature.level.1_era5_soilmoist_temp_snow<=0] <- NA # some reason NAs are shown with zero AND -247....
  Soil.temperature.level.1_era5_soilmoist_temp_snow
  summary(d$Soil.temperature.level.1_era5_soilmoist_temp_snow)
  Soil.temperature.level.1_era5_soilmoist_temp_snow <- Soil.temperature.level.1_era5_soilmoist_temp_snow/100
  names(Soil.temperature.level.1_era5_soilmoist_temp_snow) <- "Soil.temperature.level.1_era5_soilmoist_temp_snow"
  
  pred_rasters <- c(pred_rasters, Soil.temperature.level.1_era5_soilmoist_temp_snow)

  #Soil.temperature.level.1_era5_soilmoist_temp_snow <- as.data.frame(Soil.temperature.level.1_era5_soilmoist_temp_snow, xy=TRUE)
  # Topsoil temp. Temperature measured in kelvin can be converted to degrees Celsius (Â°C) by subtracting 273.15.
  
  gcs_get_object(paste0("predictors_1km/", "lst_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_1km/", "lst_", time_alt[t], ".tif"), overwrite=TRUE)
  LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality <- rast(paste0("/home/master/predictors_1km/","lst_", time_alt[t], ".tif"))
  #plot(LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality)
  LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality
  summary(d$LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality)
  LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality <- LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality/10 # changed  
  names(LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality) <- "LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality"
  #tmean_terraclimate_sites <- as.data.frame(tmean_terraclimate_sites, xy=TRUE)
  #  Mean  air temperature C degrees
  
  pred_rasters <- c(pred_rasters, LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality)
  rm(LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality)
  gc()
  

  # rm(vpd_terraclimate_sites)
  # gc()
  #vpd_terraclimate_sites <- vpd_terraclimate_sites/100  ### NO CONVERSION NEEDED 
 #vpd_terraclimate_sites <- as.data.frame(vpd_terraclimate_sites, xy=TRUE)
  #  Vapor pressure deficit kpa, both have a scale factor of   0.01 so values are really -0.001-1.4 kPA
  
  
  
  # # raster merge
  # pred_rast_dynamic1 <- c(srad_terraclimate_sites, Barrow_CO2_conc_Barrow_CO2conc,
  #                         NDVI_whittaker_constant_monthly_mean, Soil.temperature.level.1_era5_soilmoist_temp_snow, vpd_terraclimate_sites,
  #                         LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality) 
  # 
  # names(pred_rast_dynamic1) <- c("srad_terraclimate_sites", "Barrow_CO2_conc_Barrow_CO2conc",
  #                                "NDVI_whittaker_constant_monthly_mean", "Soil.temperature.level.1_era5_soilmoist_temp_snow", "vpd_terraclimate_sites",
  #                                "LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality")
  # 
  
  
  # pred_rast_dynamic1_df <- as.data.frame(pred_rast_dynamic1, xy=TRUE)
  # 
  # pred_rast_dynamic1_na <- na.omit(pred_rast_dynamic1_df)
  # str(pred_rast_dynamic1_na)
  # 
  # rm(srad_terraclimate_sites)
  # rm(Barrow_CO2_conc_Barrow_CO2conc)
  # rm(ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled)
  # rm(LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality)
  # rm(vpd_terraclimate_sites)
  # gc()
  
  
  # remove from disk too
  file.remove(paste0("/home/master/predictors_1km/", "srad_", time_alt[t], ".tif"))
  file.remove(paste0("/home/master/predictors_1km/", "co2_", time_alt[t], ".tif"))
  file.remove(paste0("/home/master/predictors_1km/", "ndvi_", time_alt[t], ".tif"))
  file.remove(paste0("/home/master/predictors_1km/", "lst_", time_alt[t], ".tif"))
  
  
  print("done")
  
  # continue with the rest of dynamic rasters...
  gcs_get_object(paste0("predictors_1km/", "snowcover_", time[t], ".tif"), saveToDisk = paste0("predictors_1km/", "snowcover_", time[t], ".tif"), overwrite=TRUE)
  Snow.cover_era5_soilmoist_temp_snow <- rast(paste0("/home/master/predictors_1km/","snowcover_", time[t], ".tif"))
  #plot(Snow.cover_era5_soilmoist_temp_snow)
  Snow.cover_era5_soilmoist_temp_snow[Snow.cover_era5_soilmoist_temp_snow< -10000] <- NA
  Snow.cover_era5_soilmoist_temp_snow <- mask(Snow.cover_era5_soilmoist_temp_snow, Soil.temperature.level.1_era5_soilmoist_temp_snow) # this only worked
  Snow.cover_era5_soilmoist_temp_snow
  Snow.cover_era5_soilmoist_temp_snow
  summary(d$Snow.cover_era5_soilmoist_temp_snow)
  Snow.cover_era5_soilmoist_temp_snow <- Snow.cover_era5_soilmoist_temp_snow/100
  names(Snow.cover_era5_soilmoist_temp_snow) <- "Snow.cover_era5_soilmoist_temp_snow"
  
  pred_rasters <- c(pred_rasters, Snow.cover_era5_soilmoist_temp_snow)
  rm(Snow.cover_era5_soilmoist_temp_snow)
  gc()
  #Snow.cover_era5_soilmoist_temp_snow <- as.data.frame(Snow.cover_era5_soilmoist_temp_snow, xy=TRUE)
  # Snow cover %
  
  gcs_get_object(paste0("predictors_1km/", "soilmoistlevel1_", time[t], ".tif"), saveToDisk = paste0("predictors_1km/", "soilmoistlevel1_", time[t], ".tif"), overwrite=TRUE)
  Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- rast(paste0("/home/master/predictors_1km/","soilmoistlevel1_", time[t], ".tif"))
  #plot(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
  Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- mask(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow, Soil.temperature.level.1_era5_soilmoist_temp_snow) # this only worked
  Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow
  summary(d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
  Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow/100
  names(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow) <- "Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow"
  
  pred_rasters <- c(pred_rasters, Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
  rm(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
  gc()
  #Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- as.data.frame(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow, xy=TRUE)
  # volumetric water content (0-1)
  
  gcs_get_object(paste0("predictors_1km/", "Percent_Tree_Cover_mod44b", substr(time[t],1, 4), ".tif"), saveToDisk = paste0("predictors_1km/", "Percent_Tree_Cover_mod44b", substr(time[t],1, 4), ".tif"), overwrite=TRUE)
  vegcover <- rast(paste0("predictors_1km/", "Percent_Tree_Cover_mod44b", substr(time[t],1, 4), ".tif"))
  vegcover <- stack(vegcover)

  Percent_Tree_Cover_MOD44B_sites <- vegcover[[1]] %>% rast()
  Percent_NonTree_Vegetation_MOD44B_sites <- vegcover[[2]] %>% rast()
  Percent_NonVegetated_MOD44B_sites <- vegcover[[3]] %>% rast()
  
  rm(vegcover)
  gc()
  
  names(Percent_Tree_Cover_MOD44B_sites) <- "Percent_Tree_Cover_MOD44B_sites"
  
  pred_rasters <- c(pred_rasters, Percent_Tree_Cover_MOD44B_sites)
  rm(Percent_Tree_Cover_MOD44B_sites)
  gc()
  
  names(Percent_NonTree_Vegetation_MOD44B_sites) <- "Percent_NonTree_Vegetation_MOD44B_sites"
  
  pred_rasters <- c(pred_rasters, Percent_NonTree_Vegetation_MOD44B_sites)
  rm(Percent_NonTree_Vegetation_MOD44B_sites)
  gc()
  
  names(Percent_NonVegetated_MOD44B_sites) <- "Percent_NonVegetated_MOD44B_sites"
  
  pred_rasters <- c(pred_rasters, Percent_NonVegetated_MOD44B_sites)
  rm(Percent_NonVegetated_MOD44B_sites)
  gc()
  

  
  setwd("/home/master/")
  gcs_get_object(paste0("predictors_1km/", "srad_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_1km/", "srad_", time_alt[t], ".tif"), overwrite=TRUE)
  Interval <- rast(paste0("/home/master/predictors_1km/", "srad_", time_alt[t], ".tif"))
  #plot(srad_terraclimate_sites)
  # classify everything to one specifric interval
  m <- c(-1000000000000000, 1000000000000000, as.numeric(substr(time[t], start=6, stop=7)))
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  Interval <- classify(Interval, rclmat, include.lowest=TRUE)
  Interval
  summary(d$Interval)
  names(Interval) <- "Interval"

  pred_rasters <- c(pred_rasters, Interval)
  rm(Interval)
  gc()
  

  gcs_get_object(paste0("predictors_1km/", "vpd_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_1km/", "vpd_", time_alt[t], ".tif"), overwrite=TRUE)
  vpd_terraclimate_sites <- rast(paste0("/home/master/predictors_1km/","vpd_", time_alt[t], ".tif"))
  #plot(vpd_terraclimate_sites)
  vpd_terraclimate_sites
  summary(d$vpd_terraclimate_sites)
  names(vpd_terraclimate_sites) <- "vpd_terraclimate_sites"
  
  pred_rasters <- c(pred_rasters, vpd_terraclimate_sites) ### TÄSSÄ NYT JOKU VIKA; TULEE TYHJÄNÄ JOS AJAA TÄMÄN AIEMMIN
  # rm(vpd_terraclimate_sites); gc() ### NÄiTÄ EI VOI POISTAA
  # file.remove(paste0("/home/master/predictors_1km/", "vpd_", time_alt[t], ".tif"))
  rm(Soil.temperature.level.1_era5_soilmoist_temp_snow)
  gc()
  
  
  # crop rasters so that there is not so much empty space - nope the tree cover is missing
  pred_rasters <- crop(pred_rasters, ext(-5e+06, 5e+06, -4.3e+06, 4.2e+06))
  
  ## now also make an additional fake raster to merge the full pred raster with that has all the model training data
  fake_rast <- crop(pred_rasters, ext(0.5e+06, 0.5e+06+(500*999.97108203905224855), -3.5e+06, -3.5e+06+(500*999.97108203905224855)))
  ext(fake_rast) <- c(0.5e+06, 0.5e+06+(500*999.97108203905224855), -4.2e+06, -4.2e+06+(500*999.97108203905224855))
  #fake_rast <- resample(fake_rast, pred_rast_static, method="near")
  # fake_rast[[1]] <- d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged # change esa cci veg type values
  fake_rast[[16]] <- d$Interval # change interval values
  gc()
  
  pred_rasters <- terra::merge(pred_rasters, fake_rast) # can be heavy to compute see https://github.com/rspatial/terra/issues/210
  rm(fake_rast)
  gc()

  pred_rasters[[16]] <- round(pred_rasters[[16]], digits=0) # had to round because the fake raster got some non-inte


  
  # rm(Interval)
  # gc()
  
  
  # # matrix approach
  # gc()
  #pred_data <- as.data.frame(pred_rasters, xy=TRUE, na.rm=TRUE)
  # pred_m <- as.matrix(pred_rasters, wide=TRUE)
  # rm(pred_rasters)
  #  # nbäiden jälkeen GiB
  # pred_df <- data.frame(pred_m)
  # pred_df_nona <- na.omit(pred_df)
  # 
  
  # print("dynamic vars pt 2 merging")
  # 
  # pred_rast_dynamic2 <- c(Snow.cover_era5_soilmoist_temp_snow, 
  #                         Percent_TreeCover_AVHRR_VCF5KYR, Percent_NonTree_Vegetation_AVHRR_VCF5KYR, Percent_NonVegetated_AVHRR_VCF5KYR, Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow, Interval) 
  # 
  # names(pred_rast_dynamic2) <- c("Snow.cover_era5_soilmoist_temp_snow", 
  #                                "Percent_TreeCover_AVHRR_VCF5KYR", "Percent_NonTree_Vegetation_AVHRR_VCF5KYR", "Percent_NonVegetated_AVHRR_VCF5KYR", "Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow", "Interval")
  # str(pred_rast_dynamic2)
  # 
  # 
  # # 
  # # pred_rast_dynamic2_df <- as.data.frame(pred_rast_dynamic2, xy=TRUE)
  # # 
  # # pred_rast_dynamic2_na <- na.omit(pred_rast_dynamic2_df)
  # # str(pred_rast_dynamic2_na)
  # 
  # rm(Snow.cover_era5_soilmoist_temp_snow)
  # rm(Percent_TreeCover_AVHRR_VCF5KYR)
  # rm(Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
  # rm(Percent_NonVegetated_AVHRR_VCF5KYR)
  # rm(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
  # 
  # 
  # gc()
  # 
  file.remove(paste0("/home/master/predictors_1km/", "vpd_", time_alt[t], ".tif"))
  file.remove(paste0("/home/master/predictors_1km/", "soilmoistlevel1_", time[t], ".tif"))
  file.remove(paste0("/home/master/predictors_1km/", "Percent_Tree_Cover_mod44b", substr(time[t],1, 4), ".tif"))
  file.remove(paste0("/home/master/predictors_1km/", "soiltemplevel1_", time[t], ".tif"))
  

  
  # # ### combine all
  # pred_rast <- c(pred_rast_static, pred_rast_dynamic1, pred_rast_dynamic2)
  
  # 
  # print("merge dynamic 1 and 2") 
  # pred_rast_dynamic <- merge(pred_rast_dynamic1_na, pred_rast_dynamic2_na, by=c("x", "y"))
  # print("merge all")
  # str(pred_rast_dynamic)
  # str(pred_rast_static_na)
  # pred_rast <- merge(pred_rast_static_na, pred_rast_dynamic, by=c("x", "y")) # rows that have NA are skipped 
  
  
  # Remove files that are not needed anymore
  #rm(pred_rast_static) # keep this in memory because we will need it later!!
  rm(pred_rast_dynamic)
  rm(pred_rast_dynamic1)
  rm(pred_rast_dynamic2)
  rm(pred_rast_dynamic1_na)
  rm(pred_rast_dynamic2_na)
  rm(pred_rast_dynamic2_df)
  rm(pred_rast_dynamic1_df)
  rm(pred_rast_static_df)
  
  gc()
  
  print("prediction data done")
  

  
  
  for (i in resp_vars) {
    
    #i <- "NEE_gC_m2"
    
    print("looping through resp vars")
    print(i)
    
    
    
    print(paste0("/home/master/abcfluxv1_modeling/results/", paste(i,  "1km_qrf", "train_loocv",sep="_"), ".rds"))
    # Load model files
    mod <- readRDS(paste0("/home/master/abcfluxv1_modeling/results/", paste(i, "1km_qrf",  "train_loocv",  sep="_"), ".rds"))
    mod2 = mod$finalModel #pull out the quantile regression object
    
    
    # # test with a smaller!!
    # pred_rasterstest <- crop(pred_rasters, ext(0.5e+06, 2e+06+(500*999.97108203905224855), -6.5e+06, -2e+06+(500*999.97108203905224855))) # tällä ennuste toimii!!
    # 
    # # data frame
    #pred <- predict(object=mod2, newdata=pred_data[1:1000000, 3:19], what=c(0.5), proximity = FALSE) # changing to 6000000 produces this error: Error in predict.randomForest(object, newdata = newdata, nodes = TRUE) : 
    # long vectors (argument 21) are not supported in .C
    # but note: predicting to this data frame without some of the factor classes works!
    
    # rasters
    # crop to two
    pred_rasters2 <- crop(pred_rasters, ext(-5e+06, 0, -4.6e+06, 4.2e+06))
    pred <- predict(model=mod2, object=pred_rasters2,  na.rm=TRUE, proximity = FALSE, filename="/home/master/temp/tempr.tif", overwrite=TRUE,  cores=10, cpkgs=c("randomForest", "quantregForest"))
    pred <- predict(model=mod2, object=pred_rasters2,  na.rm=TRUE, proximity = FALSE, filename="/home/master/temp/tempr.tif", overwrite=TRUE,  cores=1, cpkgs=c("randomForest", "quantregForest")) # changing to 1 core did not help? getting an error about Error in predict.randomForest(object, newdata = newdata, nodes = TRUE) : 
    # long vectors (argument 21) are not supported in .C
    
    # smaller domain
    pred_rasters2 <- crop(pred_rasters, ext(-5e+06, 0, -4.6e+06, 0))
    pred <- predict(model=mod2, object=pred_rasters2,  na.rm=TRUE, proximity = FALSE, filename="/home/master/temp/tempr.tif", overwrite=TRUE,  cores=10, cpkgs=c("randomForest", "quantregForest"))
    
    # full extent
    pred <- predict(model=mod2, object=pred_rasters,  na.rm=TRUE, proximity = FALSE, filename="/home/master/temp/tempr.tif", overwrite=TRUE,  cores=10, cpkgs=c("randomForest", "quantregForest")) #what = quantiles, does not work???#vector of confidence intervals to predict # https://stackoverflow.com/questions/24195805/issue-with-randomforest-long-vectors
    # 
    print("prediction to dataframe done")
    
    
    #loop through the quantiles
    for(q in 1:length(quantiles)){
      
      # q <- 1
      this_q = quantiles[[q]]
      
      #get outpath
      out_path = file.path("/home/master/predictions_8km/csv", this_q)
      
      #get the dimension of interest
      sub_pred = pred[, q]
      
      
      sub_pred <- sub_pred*10000 #multiply by 1000 so that can save in integers
      sub_pred <- round(sub_pred, digits = 0)
      
      
      # remove the predictions to the model training data which were just done to fix the random forest error
      sub_pred <- sub_pred[(nrow(modeldata1)+1):length(sub_pred)]
      
      
      
      # add cell coordinates and write out
      pred_matrix  <- data.matrix(data.frame(cbind(pred_rast[,1:2], sub_pred)))
      write.csv(pred_matrix, file.path(out_path, paste0(i,  "_8km_qrf_", time[t], ".csv")), row.names=FALSE)
      
      filename_in=paste0(out_path, "/", i,  "_8km_qrf_", time[t], ".csv")
      filename_out=paste0("gs://abcflux_modeling_files/predictions_8km/csv/", this_q, "/",  i,  "_8km_qrf_", time[t], ".csv")
      
      
      system(paste("gsutil cp", filename_in, filename_out, sep=" ")) # unix commands. use gsutils
      file.remove(filename_in)
      
    } # quantile loop
    
    print(i); " predictions done and saved"
      
    
  }  # resp var loop
  
  
} # time loop




