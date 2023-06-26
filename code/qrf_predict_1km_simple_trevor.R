


### Prepping

.libPaths("/home/master/R/x86_64-pc-linux-gnu-library/4.2")


library("dplyr")
library("caret")
library("parallel")
library("doParallel")
library("randomForest")
library("quantregForest")  
library("groupdata2")  
library(foreach)
library("sp")
library("ggplot2")
library("caret")
library("dplyr")
library("purrr")
library("raster")
library("terra")
library(stringr)
library(googleCloudStorageR)


# Terra settings
terraOptions(memfrac=0.9, tempdir = "/home/master/temp/") 
options(digits=20) 



### Load the model training data
d <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg_allsites_1km.csv")
d <- d %>% mutate_if(is.integer, as.numeric)

# Reverse gpp first so that GPP and Reco have both positive values
d$GPP_gC_m2 <- -d$GPP_gC_m2 


### Model parameters

# Predictors
Baseline_vars_1km <- c("Interval", "srad_terraclimate_sites", "vpd_terraclimate_sites",
                       
                       "Barrow_CO2_conc_Barrow_CO2conc", # atmos CO2 conc
                       
                       "ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged",
                       
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
resp_vars <- c("NEE_gC_m2", "GPP_gC_m2", "Reco_gC_m2")  




### Time periods for the monthly predictions
# Loop through average monthly rasters
time <- seq(as.Date("2001/01/01"), as.Date("2020/12/31"), "months")
time <- substr(time, 1, 7)
time <- sub("-", "_", sub("_", "", time, fixed=TRUE), fixed=TRUE)
time_alt <- gsub("_0", "_", time)

## TEMPORARY
time <- time[35:240]
time_alt <- time_alt[35:240]

# 
# ## TEMPORARY
# time <- time[7]
# time_alt <- time_alt[7]




### Load static rasters 

### Load static vars (only once)
setwd("/home/master/")


gcs_get_object("masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_northpolelambert1km_tundraboreal_simplevers.tif", saveToDisk = "masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_northpolelambert1km_tundraboreal_simplevers.tif",overwrite=TRUE)
ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- rast("masking_summary_rasters/ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_northpolelambert1km_tundraboreal_simplevers.tif")


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




print("static vars loaded")



# loop through the time periods and load dynamic data rasters - for now just extracting one

t <- 7 
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


file.remove(paste0("/home/master/predictors_1km/", "vpd_", time_alt[t], ".tif"))
file.remove(paste0("/home/master/predictors_1km/", "soilmoistlevel1_", time[t], ".tif"))
file.remove(paste0("/home/master/predictors_1km/", "Percent_Tree_Cover_mod44b", substr(time[t],1, 4), ".tif"))
file.remove(paste0("/home/master/predictors_1km/", "soiltemplevel1_", time[t], ".tif"))



pred_rast <- as.data.frame(pred_rasters, xy=TRUE, na.rm=TRUE)
rm(pred_rasters)
gc()
print("prediction data done")



### predict
i <- "NEE_gC_m2"

print("looping through resp vars")
print(i)

# test with rf
print(paste0("/home/master/abcfluxv1_modeling/results/", paste(i,  "1km_rf", "train_loocv",sep="_"), ".rds"))
# Load model files
mod <- readRDS(paste0("/home/master/abcfluxv1_modeling/results/", paste(i, "1km_rf",  "train_loocv",  sep="_"), ".rds"))
mod2 = mod$finalModel #pull out the quantile regression object


# load the model training data: we'll use this to finalize model prediction data
modeldata2 <- d[,c("Study_ID_Short", "id", i, Baseline_vars_1km)]
modeldata1 <- na.omit(modeldata2)

# editing the prediction raster based on this comment https://stackoverflow.com/questions/24829674/r-random-forest-error-type-of-predictors-in-new-data-do-not-match
# because of an error in the random forest prediction: random forest needs all the factor levels used in model training also in the prediction data

pred_rast_final <- pred_rast[, 3:ncol(pred_rast)]
pred_rast_final<-pred_rast_final[names(modeldata1)[4:20]]
pred_rast_final <- rbind(modeldata1[1:nrow(modeldata1),4:20] , pred_rast_final)

# Simplifying the land cover classes (for simplicity and because some classes have very limited amount of data)
pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==1, 31, pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)

pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==90, 30, pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)

pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==41, 160, pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)



# # Define the factors - nope, we need a numeric variable for merging the rasters later
# d$Study_ID_Short <- factor(d$Study_ID_Short)
pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- as.factor(pred_rast_final$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
pred_rast_final$Interval <- as.factor(pred_rast_final$Interval)

pred_rast_final <- na.omit(pred_rast_final)

gc()

# # tried out this parallel but it did not work
# library(itertools)
# # test parallel
# predictions <-
#   foreach(d=isplitRows(pred_rast_final, chunks=no_cores),
#           .combine=c, .packages=c("randomForest")) %dopar% {
#             predict(mod, newdata=d, proximity = FALSE)
#           }




### Set up clusters for parallel processing
#  https://stackoverflow.com/questions/44774516/parallel-processing-in-r-in-caret
# (loading the packages here, just in case, although they should already be loaded)
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
library(doParallel)
# create the cluster for caret to use
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/master/R/x86_64-pc-linux-gnu-library/4.2"))

library("caret")
library("randomForest")
library("foreach") # caret uses this
pred <- predict(object=mod, newdata=pred_rast_final, proximity = FALSE)  # use mod instead of mod2 because mod is the original caret version so that the parallel setting from the caret model will be considered

# # also tried this:
#pred <- predict.train(object=mod, newdata=pred_rast_final, proximity = FALSE)  

library(itertools)
# test parallel
predictions <-
  foreach(d=isplitRows(pred_rast_final, chunks=no_cores),
          .combine=c, .packages=c("randomForest")) %dopar% {
            predict(mod2, newdata=d, proximity = FALSE) }



sub_pred <- unname(pred)*10000 #multiply by 1000 so that can save in integers
sub_pred <- round(sub_pred, digits = 0)


# remove the predictions to the model training data which were just done to fix the random forest error
sub_pred <- sub_pred[(nrow(modeldata1)+1):length(sub_pred)]



# add cell coordinates and write out
out_path <- "/home/master/predictions_1km/"
pred_matrix  <- data.matrix(data.frame(cbind(pred_rast[,1:2], sub_pred)))
write.csv(pred_matrix, file.path(out_path, paste0(i,  "_1km_rf_", time[t], ".csv")), row.names=FALSE)



filename_in=paste0(out_path,  i,  "_1km_rf_", time[t], ".csv")
filename_out=paste0("gs://abcflux_modeling_files/predictions_1km/csv/",   i,  "_1km_rf_", time[t], ".csv")


system(paste("gsutil cp", filename_in, filename_out, sep=" ")) # unix commands. use gsutils
file.remove(filename_in)


rm(pred); rm(sub_pred); rm(pred_rast_final); rm(pred_matrix); rm(pred_rast)
print(i); " predictions done and saved"
gc()

stopCluster(cl)
registerDoSEQ()

