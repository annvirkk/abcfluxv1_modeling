

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

  print("model tuning done")

  ### Write the model out
  saveRDS(rfe_fit2, paste0("/home/master/abcfluxv1_modeling/results/", paste(flux, "8km_qrf_train_loocv_norm", sep="_"), ".rds"))



}
