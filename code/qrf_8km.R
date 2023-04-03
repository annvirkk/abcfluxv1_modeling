

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



# Define the factors
d$Study_ID_Short <- factor(d$Study_ID_Short)
d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- as.factor(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
d$Interval <- as.factor(d$Interval)



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


gcs_get_object("predictors_8km/ph.tif", saveToDisk = "predictors_8km/ph.tif", overwrite=TRUE)
PHIHOX_M_sl1_250m_ll_SoilGrids <- rast("predictors_8km/ph.tif")
#plot(PHIHOX_M_sl1_250m_ll_SoilGrids)
PHIHOX_M_sl1_250m_ll_SoilGrids
summary(d$PHIHOX_M_sl1_250m_ll_SoilGrids)
PHIHOX_M_sl1_250m_ll_SoilGrids <- PHIHOX_M_sl1_250m_ll_SoilGrids/100
#PHIHOX_M_sl1_250m_ll_SoilGrids <- as.data.frame(PHIHOX_M_sl1_250m_ll_SoilGrids, xy=TRUE)
phmin <- min(c(d$PHIHOX_M_sl1_250m_ll_SoilGrids, values(PHIHOX_M_sl1_250m_ll_SoilGrids)), na.rm=TRUE)
phmax <- max(c(d$PHIHOX_M_sl1_250m_ll_SoilGrids, values(PHIHOX_M_sl1_250m_ll_SoilGrids)), na.rm=TRUE)


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
r <- rast(files_to_download3)/10
# hist(values(rast(files_to_download3[1])/10))
# test <- d %>% filter(Interval==1); hist(test$srad_terraclimate_sites)
sradmin <- min(c(d$srad_terraclimate_sites, values(r)), na.rm=TRUE)
sradmax <- max(c(d$srad_terraclimate_sites, values(r)), na.rm=TRUE)





# tmean
select_data(var="predictors_8km/tmean")
print(files_to_download3)
r <- rast(files_to_download3)/10
# hist(values(rast(files_to_download3[1])/10))
# test <- d %>% filter(Interval==1); hist(test$tmean_terraclimate_sites)
tmeanmin <- min(c(d$tmean_terraclimate_sites, values(r)), na.rm=TRUE)
tmeanmax <- max(c(d$tmean_terraclimate_sites, values(r)), na.rm=TRUE)


# vpd
select_data(var="predictors_8km/vpd")
print(files_to_download3)
r <- rast(files_to_download3)
# hist(values(rast(files_to_download3[1])))
# test <- d %>% filter(Interval==1); hist(test$vpd_terraclimate_sites)
vpdmin <- min(c(d$vpd_terraclimate_sites, values(r)), na.rm=TRUE)
vpdmax <- max(c(d$vpd_terraclimate_sites, values(r)), na.rm=TRUE)




select_data(var="predictors_8km/co2")
print(files_to_download3)
r <- rast(files_to_download3)/1000 
# hist(values(rast(files_to_download3[1]))/1000 )
# test <- d %>% filter(Interval==1); hist(test$Barrow_CO2_conc_Barrow_CO2conc)
co2min <- min(c(d$Barrow_CO2_conc_Barrow_CO2conc , values(r)), na.rm=TRUE)
co2max <- max(c(d$Barrow_CO2_conc_Barrow_CO2conc , values(r)), na.rm=TRUE)


select_data(var="predictors_8km/ndvi_gimms")
print(files_to_download3)
r <- rast(files_to_download3)/10000 
# hist(values(rast(files_to_download3[1]))/10000 )
# test <- d %>% filter(Interval==1); hist(test$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled)
ndvimin <- min(c(d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled , values(r)), na.rm=TRUE)
ndvimax <- max(c(d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled , values(r)), na.rm=TRUE)


select_data(var="predictors_8km/soiltemplevel1")
print(files_to_download3)
r <- rast(files_to_download3)/100 
# hist(values(rast(files_to_download3[1]))/100 )
# test <- d %>% filter(Interval==1); hist(test$Soil.temperature.level.1_era5_soilmoist_temp_snow)
soiltempmin <- min(c(d$Soil.temperature.level.1_era5_soilmoist_temp_snow , values(r)), na.rm=TRUE)
soiltempmax <- max(c(d$Soil.temperature.level.1_era5_soilmoist_temp_snow , values(r)), na.rm=TRUE)




select_data(var="predictors_8km/snowcover")
print(files_to_download3)
r <- rast(files_to_download3)/100 
# hist(values(rast(files_to_download3[1]))/100 )
# test <- d %>% filter(Interval==1); hist(test$Snow.cover_era5_soilmoist_temp_snow)
snowcmin <- min(c(d$Snow.cover_era5_soilmoist_temp_snow , values(r)), na.rm=TRUE)
snowcmax <- max(c(d$Snow.cover_era5_soilmoist_temp_snow , values(r)), na.rm=TRUE)



select_data(var="predictors_8km/soilmoistlevel1")
print(files_to_download3)
r <- rast(files_to_download3)/100 
# hist(values(rast(files_to_download3[1]))/100 )
# test <- d %>% filter(Interval==1); hist(test$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow)
soilmoistmin <- min(c(d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow , values(r)), na.rm=TRUE)
soilmoistmax <- max(c(d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow , values(r)), na.rm=TRUE)





select_data(var="predictors_8km/Percent_TreeCover_VCF5KYR")
print(files_to_download3)
r <- rast(files_to_download3)
# hist(values(rast(files_to_download3[1])))
# test <- d %>% filter(Interval==1); hist(test$Percent_TreeCover_AVHRR_VCF5KYR)
treecmin <- min(c(d$Percent_TreeCover_AVHRR_VCF5KYR , values(r)), na.rm=TRUE)
treecmax <- max(c(d$Percent_TreeCover_AVHRR_VCF5KYR , values(r)), na.rm=TRUE)


select_data(var="predictors_8km/Percent_NonTree_Vegetation_VCF5KYR")
print(files_to_download3)
r <- rast(files_to_download3)
# hist(values(rast(files_to_download3[1])))
# test <- d %>% filter(Interval==1); hist(test$Percent_NonTree_Vegetation_AVHRR_VCF5KYR)
nontreecmin <- min(c(d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR , values(r)), na.rm=TRUE)
nontreecmax <- max(c(d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR , values(r)), na.rm=TRUE)


select_data(var="predictors_8km/Percent_NonVegetated_VCF5KYR")
print(files_to_download3)
r <- rast(files_to_download3)
# hist(values(rast(files_to_download3[1])))
# test <- d %>% filter(Interval==1); hist(test$Percent_NonVegetated_AVHRR_VCF5KYR)
nonvegcmin <- min(c(d$Percent_NonVegetated_AVHRR_VCF5KYR , values(r)), na.rm=TRUE)
nonvegcmax <- max(c(d$Percent_NonVegetated_AVHRR_VCF5KYR , values(r)), na.rm=TRUE)







### Do the actual normalization  (value – min) / (max – min) 
normalize <- function(x, ...) {
  return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
}
d$NEE_gC_m2 <- normalize(d$NEE_gC_m2, na.rm=TRUE)
d$Reco_gC_m2 <- normalize(d$Reco_gC_m2, na.rm=TRUE)
d$GPP_gC_m2 <- normalize(d$GPP_gC_m2, na.rm=TRUE)


d$srad_terraclimate_sites <- (d$srad_terraclimate_sites - sradmin) / (sradmax - sradmin)
d$vpd_terraclimate_sites <- (d$vpd_terraclimate_sites - vpdmin) / (vpdmax - vpdmin)
d$tmean_terraclimate_sites <- (d$tmean_terraclimate_sites - tmeanmin) / (tmeanmax - tmeanmin)
d$Barrow_CO2_conc_Barrow_CO2conc <- (d$Barrow_CO2_conc_Barrow_CO2conc - co2min) / (co2max - co2min)
d$Snow.depth_era5_soilmoist_temp_snow <- (d$Snow.depth_era5_soilmoist_temp_snow - snowdmin) / (snowdmax - snowdmin)
d$Snow.cover_era5_soilmoist_temp_snow <- (d$Snow.cover_era5_soilmoist_temp_snow - snowcmin) / (snowcmax - snowcmin)
d$Soil.temperature.level.1_era5_soilmoist_temp_snow <- (d$Soil.temperature.level.1_era5_soilmoist_temp_snow - soiltempmin) / (soiltempmax - soiltempmin)
d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow <- (d$Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow - soilmoistmin) / (soilmoistmax - soilmoistmin)
d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled <- (d$ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled - ndvimin) / (ndvimax - ndvimin)
d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m <- (d$dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m - dtmmin) / (dtmmax - dtmmin)
d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR <- (d$Percent_NonTree_Vegetation_AVHRR_VCF5KYR - nontreecmin) / (nontreecmax - nontreecmin)
d$Percent_TreeCover_AVHRR_VCF5KYR <- (d$Percent_TreeCover_AVHRR_VCF5KYR - treecmin) / (treecmax - treecmin)
d$Percent_NonVegetated_AVHRR_VCF5KYR <- (d$Percent_NonVegetated_AVHRR_VCF5KYR - nonvegcmin) / (nonvegcmax - nonvegcmin)
d$PHIHOX_M_sl1_250m_ll_SoilGrids <- (d$PHIHOX_M_sl1_250m_ll_SoilGrids - phmin) / (phmax - phmin)
d$SoilGrids_SOC_SoilGrids_SOCstock <- (d$SoilGrids_SOC_SoilGrids_SOCstock - socmin) / (socmax - socmin)
d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH <- (d$UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH - permamin) / (permamax - permamin)







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
resp_vars <- c("NEE_gC_m2",  "GPP_gC_m2", "Reco_gC_m2") # 



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




for (flux in resp_vars) {

  print(flux)
  # flux <- "NEE_gC_m2"




  modeldata3 <- d[,c("Study_ID_Short", "id", flux, Baseline_vars_20km)]
  modeldata2 <- na.omit(modeldata3)
  sapply(modeldata2, function(x) sum(is.na(x))) # no missing data
  print("20 km data set:")
  print(nrow(modeldata2))
 

  # Create a list of row indices for cross validation splits (for leave-one-fold out)
  # leave one site out used because of the highly different number of observations in some factor data levels

  set.seed(448)
  folds <- groupdata2::fold(data= modeldata2,
                            k = length(unique(modeldata2$Study_ID_Short)), id_col = "Study_ID_Short")

  # merge
  folds_data <- subset(folds, select=c( ".folds", "id"))
  modeldata2 <- full_join(modeldata2, folds_data, by="id")

  # create a row ID
  modeldata2$samplerow <- seq(1, length(modeldata2$Study_ID_Short), by=1)

  indices2 <- list()
  indices_not2 <- list()

  for (k in unique(modeldata2$.folds)) {

    # k = 1
    subs <- subset(modeldata2, .folds==k)

    # list the row ids for the fold
    sites <- list(as.integer(subs$samplerow) )

    # list other row ids that don't belong to the k fold
    sites_not <- list(as.integer(modeldata2$samplerow[!(modeldata2$samplerow %in% unlist(sites))]))

    # check that all the rows are included
    length(sites_not[[1]]) + length(sites[[1]]) == length(modeldata2$Study_ID_Short)

    # append to list (caret wants these as lists)
    # row indices for each fold, used for model evaluation
    indices2 <- append(indices2, sites)
    # row indices for model training
    indices_not2 <- append(indices_not2, sites_not)

  }


  names(indices2) <- 1:length(indices2)
  names(indices_not2) <- 1:length(indices2)



  print("model tuning and feature selection starting")



  # Model tuning parameters
  tunecontrol2 <- trainControl(
    method = "cv",  # No need to specify this because we use index column which automatically does leave-site/group-out CV
    verboseIter = FALSE,  # A logical for printing a training log. This could be set to FALSE in the final model runs
    returnData = FALSE, # A logical for saving the data. This could be set to FALSE in the final model runs
    returnResamp = "final", # A character string indicating how much of the resampled summary metrics should be saved. We only save the final model
    # p = 0.7, # For leave-group out cross-validation: the training percentage
    summaryFunction = defaultSummary, # a function to compute performance metrics across resamples.
    selectionFunction = "best", #  chooses the tuning parameter associated with the largest (or lowest for "RMSE") performance.
    index = indices_not2, # a list with elements for each resampling iteration.  needs to be integer
    indexOut = indices2, # model evaluation
    allowParallel = TRUE, # parallel processing
    savePredictions="final") # an indicator of how much of the hold-out predictions for each resample should be saved. Values can be either "all", "final", or "none".


  ### run the QRF model 
  # I'm using the non-formula method (without ~) where dummy variables for factors are not created (because trees can handle these in their own way), but e.g. for svm I'd use the formula method, so categorical variables are transformed to dummies
  
  set.seed(448)

  rfe_fit2 = train(modeldata2[,Baseline_vars_20km], modeldata2[,flux],
                  trControl = tunecontrol2, # tuning parameters
                   method="qrf", importance=TRUE)

  print("20 km rfe and tuning done")

  ### Write the model out
  saveRDS(rfe_fit2, paste0("/home/master/abcflux_modeling/results/", paste(flux, "20km_qrf_train_loocv_norm", sep="_"), ".rds"))



}
