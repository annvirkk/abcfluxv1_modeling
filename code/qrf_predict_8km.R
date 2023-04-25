


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



### Normalize model training data so that extreme values would not get as much emphasis.

### Do the actual normalization  (value – min) / (max – min) 
normalize <- function(x, ...) {
  return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
}
d$NEE_gC_m2 <- normalize(d$NEE_gC_m2, na.rm=TRUE)
d$Reco_gC_m2 <- normalize(d$Reco_gC_m2, na.rm=TRUE)
d$GPP_gC_m2 <- normalize(d$GPP_gC_m2, na.rm=TRUE)


# For environmental data, need to use the max and min values produced in the other script to quantify min and max across training data and prediction data

df <- read.csv("/home/master/abcfluxv1_modeling/results/minmax_8km.csv")

d$srad_terraclimate_sites <- (d$srad_terraclimate_sites - df$sradmin) / (df$sradmax - df$sradmin)
d$vpd_terraclimate_sites <- (d$vpd_terraclimate_sites - df$vpdmin) / (df$vpdmax - df$vpdmin)
d$tmean_terraclimate_sites <- (d$tmean_terraclimate_sites - df$tmeanmin) / (df$tmeanmax - df$tmeanmin)
d$Barrow_CO2_conc_Barrow_CO2conc <- (d$Barrow_CO2_conc_Barrow_CO2conc - df$co2min) / (df$co2max - df$co2min)
d$Snow.cover_era5_soilmoist_temp_snow <- (d$Snow.cover_era5_soilmoist_temp_snow - df$snowcmin) / (df$snowcmax - df$snowcmin)
d$Soil.temperature.level.1_era5_soilmoist_temp_snow <- (d$Soil.temperature.level.1_era5_soilmoist_temp_snow - df$soiltempmin) / (df$soiltempmax - df$soiltempmin)
d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- (d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow - df$soilmoistmin) / (df$soilmoistmax - df$soilmoistmin)
d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled <- (d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled - df$ndvimin) / (df$ndvimax - df$ndvimin)
d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- (d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m - df$dtmmin) / (df$dtmmax - df$dtmmin)
d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- (d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR - df$nontreecmin) / (df$nontreecmax - df$nontreecmin)
d$Percent_TreeCover_AVHRR_VCF5KYR <- (d$Percent_TreeCover_AVHRR_VCF5KYR - df$treecmin) / (df$treecmax - df$treecmin)
d$Percent_NonVegetated_AVHRR_VCF5KYR <- (d$Percent_NonVegetated_AVHRR_VCF5KYR - df$nonvegcmin) / (df$nonvegcmax - df$nonvegcmin)
d$PHIHOX_M_sl1_250m_ll_SoilGrids <- (d$PHIHOX_M_sl1_250m_ll_SoilGrids - df$phmin) / (df$phmax - df$phmin)
d$SoilGrids_SOC_SoilGrids_SOCstock <- (d$SoilGrids_SOC_SoilGrids_SOCstock - df$socmin) / (df$socmax - df$socmin)
d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- (d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH - df$permamin) / (df$permamax - df$permamin)



### Model parameters

# Predictors
Baseline_vars_20km <- c("srad_terraclimate_sites", "vpd_terraclimate_sites",  "tmean_terraclimate_sites",  # tmean included because don't have LST
                        "Barrow_CO2_conc_Barrow_CO2conc",    "Interval", 
                        
                        
                        "Snow.cover_era5_soilmoist_temp_snow", 
                        "Soil.temperature.level.1_era5_soilmoist_temp_snow", "Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow",
                        
                        "ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled",
                        
                        "dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m",
                        
                        "Percent_NonTree_Vegetation_AVHRR_VCF5KYR", "Percent_TreeCover_AVHRR_VCF5KYR",
                        "Percent_NonVegetated_AVHRR_VCF5KYR",# equivalent to MOD tree cover product
                        
                        "PHIHOX_M_sl1_250m_ll_SoilGrids",  "SoilGrids_SOC_SoilGrids_SOCstock",
                        
                        "UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH",
                        
                        "ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged"
                        
                        
                        
                        
)   

# check that the columns exist
Baseline_vars_20km %in% colnames(d)



### Response variables
resp_vars <- c("NEE_gC_m2",  "GPP_gC_m2", "Reco_gC_m2") 




### Time periods for the monthly predictions
# Loop through average monthly rasters
time <- seq(as.Date("1990/01/01"), as.Date("2016/12/31"), "months")
time <- substr(time, 1, 7)
time <- sub("-", "_", sub("_", "", time, fixed=TRUE), fixed=TRUE)
time_alt <- gsub("_0", "_", time)



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


gcs_get_object("masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_aggregate_northpolelambert8km_tundraboreal_attfix.tif", saveToDisk = "masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_aggregate_northpolelambert8km_tundraboreal_attfix.tif",overwrite=TRUE)
ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- rast("/home/master/masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_aggregate_northpolelambert8km_tundraboreal_attfix.tif")
# reclassifying will be done later


gcs_get_object("predictors_8km/soc.tif", saveToDisk = "predictors_8km/soc.tif", overwrite=TRUE)
SoilGrids_SOC_SoilGrids_SOCstock <-  rast("predictors_8km/soc.tif")
#plot(SoilGrids_SOC_SoilGrids_SOCstock)
SoilGrids_SOC_SoilGrids_SOCstock
summary(d$SoilGrids_SOC_SoilGrids_SOCstock)
SoilGrids_SOC_SoilGrids_SOCstock <- SoilGrids_SOC_SoilGrids_SOCstock/100
#SoilGrids_SOC_SoilGrids_SOCstock <- as.data.frame(SoilGrids_SOC_SoilGrids_SOCstock, xy=TRUE)
# Unit tonnes per ha
SoilGrids_SOC_SoilGrids_SOCstock <- (SoilGrids_SOC_SoilGrids_SOCstock- df$socmin) / (df$socmax - df$socmin)

gcs_get_object("predictors_8km/cti.tif", saveToDisk = "predictors_8km/cti.tif", overwrite=TRUE)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- rast("predictors_8km/cti.tif")
#plot(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m
summary(d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m)
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m/100
dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- (dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m - df$dtmmin) / (df$dtmmax - df$dtmmin)

#dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- as.data.frame(dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m, xy=TRUE)
# Compound topographic index, high value= high topographic moisture



gcs_get_object("predictors_8km/ph.tif", saveToDisk = "predictors_8km/ph.tif", overwrite=TRUE)
PHIHOX_M_sl1_250m_ll_SoilGrids <- rast("predictors_8km/ph.tif")
#plot(PHIHOX_M_sl1_250m_ll_SoilGrids)
PHIHOX_M_sl1_250m_ll_SoilGrids
summary(d$PHIHOX_M_sl1_250m_ll_SoilGrids)
PHIHOX_M_sl1_250m_ll_SoilGrids <- PHIHOX_M_sl1_250m_ll_SoilGrids/100
#PHIHOX_M_sl1_250m_ll_SoilGrids <- as.data.frame(PHIHOX_M_sl1_250m_ll_SoilGrids, xy=TRUE)
PHIHOX_M_sl1_250m_ll_SoilGrids <- (PHIHOX_M_sl1_250m_ll_SoilGrids - df$phmin) / (df$phmax - df$phmin)


gcs_get_object("predictors_8km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif", saveToDisk = "predictors_8km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif", overwrite=TRUE)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- rast("predictors_8km/UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH.tif")
#plot(UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH
summary(d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH)
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH/100
UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- (UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH - df$permamin) / (df$permamax - df$permamin)
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

# convert to a data frame
pred_rast_static_df <- as.data.frame(pred_rast_static, xy=TRUE)

# remove if any of the variables have NA in the pixel
pred_rast_static_na <- na.omit(pred_rast_static_df)
str(pred_rast_static_na)



print("static vars loaded")



# loop through the time periods and load dynamic data rasters
for (t in 1:length(time)) {
  
  
  
  
  # t <- 1  
  setwd("/home/master/") 
  gcs_get_object(paste0("predictors_8km/", "srad_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_8km/", "srad_", time_alt[t], ".tif"), overwrite=TRUE)
  srad_terraclimate_sites <- rast(paste0("/home/master/predictors_8km/", "srad_", time_alt[t], ".tif"))
  #plot(srad_terraclimate_sites)
  srad_terraclimate_sites
  summary(d$srad_terraclimate_sites)
  srad_terraclimate_sites <- srad_terraclimate_sites/10  
  srad_terraclimate_sites <- (srad_terraclimate_sites - df$sradmin) / (df$sradmax - df$sradmin)  #srad_terraclimate_sites <- as.data.frame(srad_terraclimate_sites, xy=TRUE)
  # Downward surface shortwave radiation. Unit W/m2. Both need to be divided by 10 to get to the original scale
  
  
  gcs_get_object(paste0("predictors_8km/", "co2_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_8km/", "co2_", time_alt[t], ".tif"), overwrite=TRUE)
  Barrow_CO2_conc_Barrow_CO2conc <- rast(paste0("/home/master/predictors_8km/", "co2_", time_alt[t], ".tif"))
  #plot(Barrow_CO2_conc_Barrow_CO2conc)
  Barrow_CO2_conc_Barrow_CO2conc
  summary(d$Barrow_CO2_conc_Barrow_CO2conc)
  Barrow_CO2_conc_Barrow_CO2conc <- Barrow_CO2_conc_Barrow_CO2conc/1000  
  Barrow_CO2_conc_Barrow_CO2conc <- (Barrow_CO2_conc_Barrow_CO2conc - df$co2min) / (df$co2max - df$co2min)  #Barrow_CO2_conc_Barrow_CO2conc <- as.data.frame(Barrow_CO2_conc_Barrow_CO2conc, xy=TRUE)
  # atm CO2 concentrations in ppm
  
  
  gcs_get_object(paste0("predictors_8km/", "ndvi_gimms_", time[t], ".tif"), saveToDisk = paste0("predictors_8km/", "ndvi_gimms_", time[t], ".tif"), overwrite=TRUE)
  ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled <- rast(paste0("/home/master/predictors_8km/","ndvi_gimms_", time[t], ".tif")) 
  #plot(NDVI_whittaker_constant_monthly_mean)
  ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled
  summary(d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled)
  ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled <- ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled/10000 
  ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled <- (ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled  - df$ndvimin) / (df$ndvimax - df$ndvimin)
  #NDVI_whittaker_constant_monthly_mean <- as.data.frame(NDVI_whittaker_constant_monthly_mean, xy=TRUE)
  
  
  gcs_get_object(paste0("predictors_8km/", "soiltemplevel1_", time[t], ".tif"), saveToDisk = paste0("predictors_8km/", "soiltemplevel1_", time[t], ".tif"), overwrite=TRUE)
  Soil.temperature.level.1_era5_soilmoist_temp_snow <- rast(paste0("/home/master/predictors_8km/","soiltemplevel1_", time[t], ".tif"))
  #plot(Soil.temperature.level.1_era5_soilmoist_temp_snow)
  Soil.temperature.level.1_era5_soilmoist_temp_snow
  summary(d$Soil.temperature.level.1_era5_soilmoist_temp_snow)
  Soil.temperature.level.1_era5_soilmoist_temp_snow <- Soil.temperature.level.1_era5_soilmoist_temp_snow/100
  Soil.temperature.level.1_era5_soilmoist_temp_snow <- (Soil.temperature.level.1_era5_soilmoist_temp_snow - df$soiltempmin) / (df$soiltempmax - df$soiltempmin)
  #Soil.temperature.level.1_era5_soilmoist_temp_snow <- as.data.frame(Soil.temperature.level.1_era5_soilmoist_temp_snow, xy=TRUE)
  # Topsoil temp. Temperature measured in kelvin can be converted to degrees Celsius (Â°C) by subtracting 273.15.
  
  gcs_get_object(paste0("predictors_8km/", "tmean_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_8km/", "tmean_", time_alt[t], ".tif"), overwrite=TRUE)
  tmean_terraclimate_sites <- rast(paste0("/home/master/predictors_8km/","tmean_", time_alt[t], ".tif"))
  #plot(vpd_terraclimate_sites)
  tmean_terraclimate_sites
  summary(d$tmean_terraclimate_sites)
  tmean_terraclimate_sites <- tmean_terraclimate_sites/10  
  tmean_terraclimate_sites <- (tmean_terraclimate_sites - df$tmeanmin) / (df$tmeanmax - df$tmeanmin)
  #tmean_terraclimate_sites <- as.data.frame(tmean_terraclimate_sites, xy=TRUE)
  #  Mean annual air temperature C degrees
  
  gcs_get_object(paste0("predictors_8km/", "vpd_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_8km/", "vpd_", time_alt[t], ".tif"), overwrite=TRUE)
  vpd_terraclimate_sites <- rast(paste0("/home/master/predictors_8km/","vpd_", time_alt[t], ".tif"))
  #plot(vpd_terraclimate_sites)
  vpd_terraclimate_sites
  summary(d$vpd_terraclimate_sites)
  #vpd_terraclimate_sites <- vpd_terraclimate_sites/100  ### NO CONVERSION NEEDED 
  vpd_terraclimate_sites <- (vpd_terraclimate_sites- df$vpdmin) / (df$vpdmax - df$vpdmin)  #vpd_terraclimate_sites <- as.data.frame(vpd_terraclimate_sites, xy=TRUE)
  #  Vapor pressure deficit kpa, both have a scale factor of   0.01 so values are really -0.001-1.4 kPA
  
  
  
  # raster merge
  pred_rast_dynamic1 <- c(srad_terraclimate_sites, Barrow_CO2_conc_Barrow_CO2conc,
                          ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled, Soil.temperature.level.1_era5_soilmoist_temp_snow, vpd_terraclimate_sites,
                          tmean_terraclimate_sites) 
  
  names(pred_rast_dynamic1) <- c("srad_terraclimate_sites", "Barrow_CO2_conc_Barrow_CO2conc",
                                 "ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled", "Soil.temperature.level.1_era5_soilmoist_temp_snow", "vpd_terraclimate_sites",
                                 "tmean_terraclimate_sites")
  
  
  
  pred_rast_dynamic1_df <- as.data.frame(pred_rast_dynamic1, xy=TRUE)
  
  pred_rast_dynamic1_na <- na.omit(pred_rast_dynamic1_df)
  str(pred_rast_dynamic1_na)
  
  rm(srad_terraclimate_sites)
  rm(Barrow_CO2_conc_Barrow_CO2conc)
  rm(ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled)
  rm(Soil.temperature.level.1_era5_soilmoist_temp_snow)
  rm(tmean_terraclimate_sites)
  rm(vpd_terraclimate_sites)
  gc()
  
  
  # remove from disk too
  file.remove(paste0("/home/master/predictors_8km/", "srad_", time_alt[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "co2_", time_alt[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "ndvi_gimms_", time[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "soiltemplevel1_", time[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "tmean_", time_alt[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "vpd_", time_alt[t], ".tif"))
  
  
  print("done")
  
  # continue with the rest of dynamic rasters...
  gcs_get_object(paste0("predictors_8km/", "snowcover_", time[t], ".tif"), saveToDisk = paste0("predictors_8km/", "snowcover_", time[t], ".tif"), overwrite=TRUE)
  Snow.cover_era5_soilmoist_temp_snow <- rast(paste0("/home/master/predictors_8km/","snowcover_", time[t], ".tif"))
  #plot(Snow.cover_era5_soilmoist_temp_snow)
  Snow.cover_era5_soilmoist_temp_snow
  summary(d$Snow.cover_era5_soilmoist_temp_snow)
  Snow.cover_era5_soilmoist_temp_snow <- Snow.cover_era5_soilmoist_temp_snow/100
  Snow.cover_era5_soilmoist_temp_snow <- (Snow.cover_era5_soilmoist_temp_snow- df$snowcmin) / (df$snowcmax - df$snowcmin)  
  #Snow.cover_era5_soilmoist_temp_snow <- as.data.frame(Snow.cover_era5_soilmoist_temp_snow, xy=TRUE)
  # Snow cover %
  
  gcs_get_object(paste0("predictors_8km/", "soilmoistlevel1_", time[t], ".tif"), saveToDisk = paste0("predictors_8km/", "soilmoistlevel1_", time[t], ".tif"), overwrite=TRUE)
  Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- rast(paste0("/home/master/predictors_8km/","soilmoistlevel1_", time[t], ".tif"))
  #plot(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
  Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow
  summary(d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
  Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow/100
  Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- (Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow- df$soilmoistmin) / (df$soilmoistmax - df$soilmoistmin)  #Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- as.data.frame(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow, xy=TRUE)
  # volumetric water content (0-1)
  
  
  # AVHRR fields lacking data from 1994 and 2000 - use data from the previous year in this case
  if (substr(time[t],1, 4)==1994) {
    
    gcs_get_object(paste0("predictors_8km/", "Percent_TreeCover_VCF5KYR_", 1993, "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_TreeCover_VCF5KYR_", 1993, "001.tif"), overwrite=TRUE)
    Percent_TreeCover_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_TreeCover_VCF5KYR_", 1993, "001.tif"))
    print(Percent_TreeCover_AVHRR_VCF5KYR)
    #plot(Percent_TreeCover_AVHRR_VCF5KYR)
    Percent_TreeCover_AVHRR_VCF5KYR
    summary(d$Percent_TreeCover_AVHRR_VCF5KYR)
    Percent_TreeCover_AVHRR_VCF5KYR <- (Percent_TreeCover_AVHRR_VCF5KYR- df$treecmin) / (df$treecmax - df$treecmin)    #Percent_TreeCover_AVHRR_VCF5KYR <- as.data.frame(Percent_TreeCover_AVHRR_VCF5KYR, xy=TRUE)
    
    gcs_get_object(paste0("predictors_8km/", "Percent_NonTree_Vegetation_VCF5KYR_", 1993, "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_NonTree_Vegetation_VCF5KYR_", 1993, "001.tif"), overwrite=TRUE)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_NonTree_Vegetation_VCF5KYR_", 1993, "001.tif"))
    print(Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    #plot(Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR
    summary(d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- (Percent_NonTree_Vegetation_AVHRR_VCF5KYR- df$nontreecmin) / (df$nontreecmax - df$nontreecmin)    #Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- as.data.frame(Percent_NonTree_Vegetation_AVHRR_VCF5KYR, xy=TRUE)
    
    gcs_get_object(paste0("predictors_8km/", "Percent_NonVegetated_VCF5KYR_", 1993, "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_NonVegetated_VCF5KYR_", 1993, "001.tif"), overwrite=TRUE)
    Percent_NonVegetated_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_NonVegetated_VCF5KYR_", 1993, "001.tif"))
    print(Percent_NonVegetated_AVHRR_VCF5KYR)
    #plot(Percent_NonVegetated_AVHRR_VCF5KYR)
    Percent_NonVegetated_AVHRR_VCF5KYR
    summary(d$Percent_NonVegetated_AVHRR_VCF5KYR)
    Percent_NonVegetated_AVHRR_VCF5KYR <- (Percent_NonVegetated_AVHRR_VCF5KYR- df$nonvegcmin) / (df$nonvegcmax - df$nonvegcmin)    #Percent_NonVegetated_AVHRR_VCF5KYR <- as.data.frame(Percent_NonVegetated_AVHRR_VCF5KYR, xy=TRUE)
    
    
  } else if (substr(time[t],1, 4)==2000) {
    gcs_get_object(paste0("predictors_8km/", "Percent_TreeCover_VCF5KYR_", 1999, "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_TreeCover_VCF5KYR_", 1999, "001.tif"), overwrite=TRUE)
    Percent_TreeCover_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_TreeCover_VCF5KYR_", 1999, "001.tif"))
    print(Percent_TreeCover_AVHRR_VCF5KYR)
    #plot(Percent_TreeCover_AVHRR_VCF5KYR)
    Percent_TreeCover_AVHRR_VCF5KYR
    summary(d$Percent_TreeCover_AVHRR_VCF5KYR)
    Percent_TreeCover_AVHRR_VCF5KYR <- (Percent_TreeCover_AVHRR_VCF5KYR- df$treecmin) / (df$treecmax - df$treecmin)    #Percent_TreeCover_AVHRR_VCF5KYR <- as.data.frame(Percent_TreeCover_AVHRR_VCF5KYR, xy=TRUE)
    
    gcs_get_object(paste0("predictors_8km/", "Percent_NonTree_Vegetation_VCF5KYR_", 1999, "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_NonTree_Vegetation_VCF5KYR_", 1999, "001.tif"), overwrite=TRUE)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_NonTree_Vegetation_VCF5KYR_", 1999, "001.tif"))
    print(Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    #plot(Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR
    summary(d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- (Percent_NonTree_Vegetation_AVHRR_VCF5KYR- df$nontreecmin) / (df$nontreecmax - df$nontreecmin)
    #Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- as.data.frame(Percent_NonTree_Vegetation_AVHRR_VCF5KYR, xy=TRUE)
    
    gcs_get_object(paste0("predictors_8km/", "Percent_NonVegetated_VCF5KYR_", 1999, "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_NonVegetated_VCF5KYR_", 1999, "001.tif"), overwrite=TRUE)
    Percent_NonVegetated_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_NonVegetated_VCF5KYR_", 1999, "001.tif"))
    print(Percent_NonVegetated_AVHRR_VCF5KYR)
    #plot(Percent_NonVegetated_AVHRR_VCF5KYR)
    Percent_NonVegetated_AVHRR_VCF5KYR
    summary(d$Percent_NonVegetated_AVHRR_VCF5KYR)
    Percent_NonVegetated_AVHRR_VCF5KYR <- (Percent_NonVegetated_AVHRR_VCF5KYR- df$nonvegcmin) / (df$nonvegcmax - df$nonvegcmin)    #Percent_NonVegetated_AVHRR_VCF5KYR <- as.data.frame(Percent_NonVegetated_AVHRR_VCF5KYR, xy=TRUE)
    
    
  } else {
    
    gcs_get_object(paste0("predictors_8km/", "Percent_TreeCover_VCF5KYR_", substr(time[t],1, 4), "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_TreeCover_VCF5KYR_", substr(time[t],1, 4), "001.tif"), overwrite=TRUE)
    Percent_TreeCover_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_TreeCover_VCF5KYR_", substr(time[t],1, 4), "001.tif"))
    print(Percent_TreeCover_AVHRR_VCF5KYR)
    #plot(Percent_TreeCover_AVHRR_VCF5KYR)
    Percent_TreeCover_AVHRR_VCF5KYR
    summary(d$Percent_TreeCover_AVHRR_VCF5KYR)
    Percent_TreeCover_AVHRR_VCF5KYR <- (Percent_TreeCover_AVHRR_VCF5KYR- df$treecmin) / (df$treecmax - df$treecmin)    #Percent_TreeCover_AVHRR_VCF5KYR <- as.data.frame(Percent_TreeCover_AVHRR_VCF5KYR, xy=TRUE)
    
    gcs_get_object(paste0("predictors_8km/", "Percent_NonTree_Vegetation_VCF5KYR_", substr(time[t],1, 4), "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_NonTree_Vegetation_VCF5KYR_", substr(time[t],1, 4), "001.tif"), overwrite=TRUE)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_NonTree_Vegetation_VCF5KYR_", substr(time[t],1, 4), "001.tif"))
    print(Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    #plot(Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR
    summary(d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
    Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- (Percent_NonTree_Vegetation_AVHRR_VCF5KYR- df$nontreecmin) / (df$nontreecmax - df$nontreecmin)    #Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- as.data.frame(Percent_NonTree_Vegetation_AVHRR_VCF5KYR, xy=TRUE)
    
    gcs_get_object(paste0("predictors_8km/", "Percent_NonVegetated_VCF5KYR_", substr(time[t],1, 4), "001.tif"), saveToDisk = paste0("predictors_8km/", "Percent_NonVegetated_VCF5KYR_", substr(time[t],1, 4), "001.tif"), overwrite=TRUE)
    Percent_NonVegetated_AVHRR_VCF5KYR <- rast(paste0("/home/master/predictors_8km/","Percent_NonVegetated_VCF5KYR_", substr(time[t],1, 4), "001.tif"))
    print(Percent_NonVegetated_AVHRR_VCF5KYR)
    #plot(Percent_NonVegetated_AVHRR_VCF5KYR)
    Percent_NonVegetated_AVHRR_VCF5KYR
    summary(d$Percent_NonVegetated_AVHRR_VCF5KYR)
    Percent_NonVegetated_AVHRR_VCF5KYR <- (Percent_NonVegetated_AVHRR_VCF5KYR- df$nonvegcmin) / (df$nonvegcmax - df$nonvegcmin)    #Percent_NonVegetated_AVHRR_VCF5KYR <- as.data.frame(Percent_NonVegetated_AVHRR_VCF5KYR, xy=TRUE)
    
    
    
  }
  
  
  setwd("/home/master/")
  gcs_get_object(paste0("predictors_8km/", "srad_", time_alt[t], ".tif"), saveToDisk = paste0("predictors_8km/", "srad_", time_alt[t], ".tif"), overwrite=TRUE)
  Interval <- rast(paste0("/home/master/predictors_8km/", "srad_", time_alt[t], ".tif"))
  #plot(srad_terraclimate_sites)
  # classify everything to one specifric interval
  m <- c(-1000000000000000, 1000000000000000, as.numeric(substr(time[t], start=6, stop=7)))
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  Interval <- classify(Interval, rclmat, include.lowest=TRUE)
  Interval
  summary(d$Interval)
  
  
  
  
  
  
  
  
  gc()
  
  
  
  
  print("dynamic vars pt 2 merging")
  
  pred_rast_dynamic2 <- c(Snow.cover_era5_soilmoist_temp_snow, 
                          Percent_TreeCover_AVHRR_VCF5KYR, Percent_NonTree_Vegetation_AVHRR_VCF5KYR, Percent_NonVegetated_AVHRR_VCF5KYR, Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow, Interval) 
  
  names(pred_rast_dynamic2) <- c("Snow.cover_era5_soilmoist_temp_snow", 
                                 "Percent_TreeCover_AVHRR_VCF5KYR", "Percent_NonTree_Vegetation_AVHRR_VCF5KYR", "Percent_NonVegetated_AVHRR_VCF5KYR", "Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow", "Interval")
  str(pred_rast_dynamic2)
  
  
  
  pred_rast_dynamic2_df <- as.data.frame(pred_rast_dynamic2, xy=TRUE)
  
  pred_rast_dynamic2_na <- na.omit(pred_rast_dynamic2_df)
  str(pred_rast_dynamic2_na)
  
  rm(Snow.cover_era5_soilmoist_temp_snow)
  rm(Percent_TreeCover_AVHRR_VCF5KYR)
  rm(Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
  rm(Percent_NonVegetated_AVHRR_VCF5KYR)
  rm(Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
  
  gc()
  
  
  file.remove(paste0("/home/master/predictors_8km/", "vpd_", time_alt[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "soilmoistlevel1_", time[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "Percent_TreeCover_VCF5KYR_", time[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "Percent_NonTree_Vegetation_VCF5KYR_", time[t], ".tif"))
  file.remove(paste0("/home/master/predictors_8km/", "Percent_NonVegetated_VCF5KYR_", time[t], ".tif"))

  
  ### combine all
  
  print("merge dynamic 1 and 2") 
  pred_rast_dynamic <- merge(pred_rast_dynamic1_na, pred_rast_dynamic2_na, by=c("x", "y"))
  print("merge all")
  str(pred_rast_dynamic)
  str(pred_rast_static_na)
  pred_rast <- merge(pred_rast_static_na, pred_rast_dynamic, by=c("x", "y")) # rows that have NA are skipped 
  
  
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
    
    
    # load the model training data: we'll use this to finalize model prediction data
    modeldata2 <- d[,c("Study_ID_Short", "id", i, Baseline_vars_20km)]
    modeldata1 <- na.omit(modeldata2)

    # editing the prediction raster based on this comment https://stackoverflow.com/questions/24829674/r-random-forest-error-type-of-predictors-in-new-data-do-not-match
    # because of an error in the random forest prediction: random forest needs all the factor levels used in model training also in the prediction data
    
    pred_rast_final <- pred_rast[, 3:ncol(pred_rast)]
    pred_rast_final<-pred_rast_final[names(modeldata1)[4:20]]
    pred_rast_final <- rbind(modeldata1[1:nrow(modeldata1),4:20] , pred_rast_final)
    
    # # then convert to factors, this is a must
    pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==1, 31, pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
    
    pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==90, 30, pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
    
    pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==41, 160, pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
    
    pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- as.factor(pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
    
    pred_rast_final$Interval <- as.factor(pred_rast_final$Interval)
    
    
    
    print(paste0("/home/master/abcfluxv1_modeling/results/", paste(i,  "8km_qrf", "train_loocv_norm",sep="_"), ".rds"))
    # Load model files
    mod <- readRDS(paste0("/home/master/abcfluxv1_modeling/results/", paste(i, "8km_qrf",  "train_loocv_norm",  sep="_"), ".rds"))
    mod2 = mod$finalModel #pull out the quantile regression object
    
    
    pred <- predict(mod2, newdata=pred_rast_final, what = quantiles) #vector of confidence intervals to predict
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




