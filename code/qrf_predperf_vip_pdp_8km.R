

.libPaths("/home/master/R/x86_64-pc-linux-gnu-library/4.2")

library("caret")
library("vip")
library("pdp")
library(Metrics)
library("tdr") # metrica would be even better https://cran.r-project.org/web/packages/metrica/vignettes/available_metrics_regression.html
library("ggplot2")
library("viridis")
library("tidyr")
library("dplyr")
library("ggridges")
library(googleCloudStorageR)


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








# merge some site-level data back
d_orig <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg.csv")
d_orig <- subset(d_orig, select=c("Study_ID_Short", "Biome", "Disturbance", "Country")) %>% unique()
duplicated(d_orig$Study_ID_Short)

# test <- merge(d, d_orig, by="Study_ID_Short")
# unique(d$Study_ID_Short) %in% unique(test$Study_ID_Short)
# unique(test$Study_ID_Short)

d <- merge(d, d_orig, by="Study_ID_Short")









### Normalize model training data so that extreme values would not get as much emphasis.

### Do the actual normalization  (value – min) / (max – min) 
normalize <- function(x, ...) {
  return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
}
d$NEE_gC_m2_scaled <- normalize(d$NEE_gC_m2, na.rm=TRUE)



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




d$Study_ID_Short[d$Study_ID_Short=="L\xf3pez-Blanco_GL-NuF_tower1"] <- "Lopez-Blanco_GL-NuF_tower1"
d$Study_ID_Short[d$Study_ID_Short=="L\xf3pez-Blanco_GL-ZaF_tower1"] <- "Lopez-Blanco_GL-Zaf_tower1"

# Define the factors
d$Study_ID_Short <- factor(d$Study_ID_Short)
d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- as.factor(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
d$Interval <- as.factor(d$Interval)



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





### Set folder for results
setwd("/home/master/abcfluxv1_modeling/figures")



# For figures

theme_pub <- theme_bw() + theme(panel.border=element_rect(size=1, colour="black"),
                                axis.text=element_text(size=14, face="bold"),
                                plot.title = element_text(size = 14, face = "bold"),
                                axis.title=element_text(size=14, face="bold"), 
                                plot.subtitle=element_text(size=14, face="bold", color="black"), 
                                strip.text.x = element_text(size = 14, face="bold"),
                                legend.text=element_text(size=14, face="bold"), legend.title=element_text(size=14))



for (i in resp_vars) {
  
  #i <- "NEE_gC_m2"
  #i <- "GPP_gC_m2"
  #i <- "Reco_gC_m2"
  
  # Add data

  modeldata2 <- d[,c("Study_ID_Short", "id", i, Baseline_vars_20km)]
  modeldata2 <- na.omit(modeldata2) 

  # create a row ID
  modeldata2$samplerow <- seq(1, length(modeldata2$Study_ID_Short), by=1)
  
  
  # merge back other information
  modeldata22 <- merge(modeldata2[ , !(names(modeldata2) %in%  c("Study_ID_Short",  i, Baseline_vars_20km))], d, by="id")
  
  
  
    # i <- "NEE_gC_m2"
    # i <- "GPP_gC_m2"
    # i <- "Reco_gC_m2" 
    
    

      
  print(paste0("/home/master/abcfluxv1_modeling/results/", paste(i,  "8km_qrf", "train_loocv_norm",sep="_"), ".rds"))
  # Load model files
  mod <- readRDS(paste0("/home/master/abcfluxv1_modeling/results/", paste(i, "8km_qrf",  "train_loocv_norm",  sep="_"), ".rds"))
  mod2 = mod$finalModel #pull out the quantile regression object
      

  ### Predictive performance plots ###
  
  
  # Define x and y lab titles for the plot
  
  
  if (i=="GPP_gC_m2") {
    ylab = expression(paste("Observed GPP g C m"^{-2}, month^{-1}))
    xlab = expression(paste("Predicted GPP g C m"^{-2}, month^{-1}, " (8 km model)"))
    
  }
  
  if (i=="NEE_gC_m2") {
    ylab = expression(paste("Observed NEE g C m"^{-2}, month^{-1}))
    xlab = expression(paste("Predicted NEE g C m"^{-2}, month^{-1}, " (8 km model)"))
    
  }
  
  if (i=="Reco_gC_m2") {
    ylab = expression(paste("Observed Reco g C m"^{-2}, month^{-1}))
    xlab = expression(paste("Predicted Reco g C m"^{-2}, month^{-1}, " (8 km model)"))
    
  }
  
  
  
  # # Extract the final model details
  preds <- mod$pred %>%
    data.frame()
  
  
  
  
  # Merge

  obspred <- merge(modeldata22, preds, by.x="samplerow", by.y="rowIndex")
  # plot(obspred$NEE_gC_m2_scaled, obspred$obs) # yes, identical
  # cor.test(obspred$NEE_gC_m2, obspred$obs)
  # 
  
  
  # Rescale
  
  minval <- ifelse(i=="NEE_gC_m2", min(d$NEE_gC_m2, na.rm=TRUE), NA)
  minval <- ifelse(i=="GPP_gC_m2", min(d$GPP_gC_m2, na.rm=TRUE), minval)
  minval <- ifelse(i=="Reco_gC_m2", min(d$Reco_gC_m2, na.rm=TRUE), minval)
  
  
  maxval <- ifelse(i=="NEE_gC_m2", max(d$NEE_gC_m2, na.rm=TRUE), NA)
  maxval <- ifelse(i=="GPP_gC_m2", max(d$GPP_gC_m2, na.rm=TRUE), maxval)
  maxval <- ifelse(i=="Reco_gC_m2", max(d$Reco_gC_m2, na.rm=TRUE), maxval)
  
  # normalization values back
  # every third column, do this

  obspred[, c("pred", "obs")] <- obspred[, c("pred", "obs")]  /1 * (maxval - minval) + minval
  

    

  
  ggplot(obspred) + geom_point(aes(x=Interval, y=obs))  + geom_point(aes(x=Interval, y=pred), col="red") + facet_wrap(~Biome) + theme_pub
  
  ggplot(obspred) + geom_point(aes(x=Interval, y=obs, col=Disturbance))  + facet_wrap(~Biome) + theme_pub
  
  # First plot scatterplots for each variable based on individual models
  # Max and min of several columns 
  scale_max <- max(c(obspred$obs, obspred$pred))
  scale_min <- min(c(obspred$obs, obspred$pred))
  
  
  # Merge pred perf
  # not sure how to access the results from the best model! This is a shortcut: mod$results[which.min(mod$results[, "RMSE"]), ]
  # note that it looks like predictive performances are calculated separately for each fold (group) after which a mean is calculated. See: mean(mod$resample$Rsquared, na.rm=TRUE)
  # r_stats <- data.frame(cbind(RMSE=c(mod$results[which.min(mod$results[, "RMSE"]), ]$RMSE),
  #                             Rsquared=(mod$results[which.min(mod$results[, "RMSE"]), ]$Rsquared),
  #                             MAE=c(mod$results[which.min(mod$results[, "RMSE"]), ]$MAE)))
  # r_stats$RMSE <- as.character(r_stats$RMSE) %>% as.numeric()
  # r_stats$Rsquared <- as.character(r_stats$Rsquared) %>% as.numeric()
  # r_stats$MAE <- as.character(r_stats$MAE) %>% as.numeric()
  
  # need to calculate these by myself to control for re-scaling
  Rsquared <- cor(obspred$obs, obspred$pred)^2
  rmse <- rmse(obspred$obs, obspred$pred)
  mae <- mae(obspred$obs, obspred$pred)
  mbe <- tdStats(obspred$pred,obspred$obs,  functions = "mbe") # Negative values indicate overestimation
  r_stats <- data.frame(Rsquared=Rsquared, RMSE=rmse, MAE=mae, MBE=mbe)
  
  # how many sites?
  
  length(unique(obspred$Study_ID_Short))
  
  ### Plot XGBOOST based on categorical information that we have
  
  # Colored by biome
  p1 <- ggplot(obspred, aes(x=pred, y=obs, colour=factor(Biome))) +
    geom_point(shape = 16,  size=3) + geom_abline(slope = 1) + 
    
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    
    theme_pub  + theme(legend.title=element_blank()) + scale_colour_viridis(discrete=TRUE) + xlab(xlab) + ylab(ylab)+ 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
  
  # Print out
  setwd("/home/master/abcfluxv1_modeling/figures/")
  print(p1)
  dev.copy(png, paste(i, "qrf_8km_train_loocv_predperf_biome.png", sep="_"), width=500, height=400)
  
  dev.off()
  
  
  # Colored by veg type
  p2 <- ggplot(obspred, aes(x=pred, y=obs, colour=factor(ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged))) +
    geom_point(shape = 16,  size=3) + geom_abline(slope = 1) + 
    
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    
    theme_pub  + theme(legend.title=element_blank()) + scale_colour_viridis(discrete=TRUE) + xlab(xlab) + ylab(ylab)+  
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
  
  print(p2)
  dev.copy(png, paste(i, "qrf_8km_train_loocv_predperf_vegtype.png", sep="_"), width=500, height=400)
  dev.off()
  
  
  # Colored by disturbance
  obspred$Disturbance <- as.character(obspred$Disturbance)
  obspred$Disturbance <- ifelse(is.na(obspred$Disturbance), "No", obspred$Disturbance)
  obspred$Disturbance <- ifelse(obspred$Disturbance=="No", "NA/No", obspred$Disturbance)
  obspred$Disturbance <- factor(obspred$Disturbance, levels=c( "NA/No", "Thermokarst", "Drainage",  "Fire", "Harvest", "Larval Outbreak"))
  
  p3 <- ggplot(obspred, aes(x=pred, y=obs, colour=Disturbance)) +
    geom_point(shape = 16,  size=3, alpha=0.7) + geom_abline(slope = 1) + 
    
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    
    theme_pub  + theme(legend.title=element_blank())  + scale_colour_viridis(discrete=TRUE) + xlab(xlab) + ylab(ylab)+ 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
  
  print(p3)
  dev.copy(png, paste(i, "qrf_8km_train_loocv_predperf_disturbance.png", sep="_"), width=600, height=400)
  dev.off()
  
  
  
  # Colored by flux method
  p5 <- ggplot(obspred, aes(x=pred, y=obs, colour=factor(Flux_method))) +
    geom_point(shape = 16,  size=3) + geom_abline(slope = 1) + 
    
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    
    theme_pub  + theme(legend.title=element_blank())  + scale_colour_viridis(discrete=TRUE) + xlab(xlab) + ylab(ylab)+  
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
  
  print(p5)
  dev.copy(png, paste(i, "qrf_8km_train_loocv_predperf_fluxmethod.png", sep="_"), width=500, height=400)
  dev.off()
  
  
  # Colored by flux method detail
  p6 <- ggplot(obspred, aes(x=pred, y=obs, colour=factor(Flux_method_detail))) +
    geom_point(shape = 16,  size=3) + geom_abline(slope = 1) + 
    
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    
    theme_pub  + theme(legend.title=element_blank()) + scale_colour_viridis(discrete=TRUE) + xlab(xlab) + ylab(ylab)+ 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
  
  
  print(p6)
  dev.copy(png, paste(i, "qrf_8km_train_loocv_predperf_fluxmethoddetail.png", sep="_"), width=500, height=400)
  dev.off()
  
  
  
  
  
  # Color only the sites with outlier observations
  obspred$Study_ID_figure <- ifelse(obspred$obs-obspred$pred>80 | obspred$pred-obspred$obs>80, as.character(obspred$Study_ID_Short), "ok")

  # obspred$Study_ID_figure <- ifelse(obspred$obs-obspred$pred< -0.3 | obspred$pred-obspred$obs>0.75, as.character(obspred$Study_ID_Short), "ok")
  # obspred$Study_ID_figure <- ifelse(obspred$obs-obspred$pred< -0.2 | obspred$pred-obspred$obs>0.2, as.character(obspred$Study_ID_Short), "ok")
  p8 <- ggplot(obspred, aes(x=pred, y=obs, colour=factor(Study_ID_figure))) +
    geom_point(shape = 16,  size=3) + geom_abline(slope = 1) + 
    
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    
    theme_pub  + theme(legend.title=element_blank())  + scale_colour_viridis(discrete=TRUE) + xlab(xlab) + ylab(ylab)+ 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
  
  print(p8)
  dev.copy(png, paste(i, "qrf_8km_train_loocv_predperf_outliers.png", sep="_"), width=500, height=400)
  dev.off()
  
  
  
  # Color only the sites with months - asked by Brendan
  p9 <- ggplot(obspred, aes(x=pred, y=obs, colour=factor(Interval))) +
    geom_point(shape = 16,  size=3, alpha=0.5) + geom_abline(slope = 1) + 
    
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    
    theme_pub  + theme(legend.title=element_blank())  + scale_colour_viridis(discrete=TRUE) + xlab(xlab) + ylab(ylab)+  
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
  
  print(p9)
  dev.copy(png, paste(i, "8km_train_loocv_predperf_months.png", sep="_"), width=500, height=400)
  dev.off()
  
  
  
  p9 <- ggplot(obspred, aes(x=pred, y=obs, colour=Interval)) +
    geom_point(shape = 16,  size=3) + geom_abline(slope = 1) + 
    
    theme_pub  + theme(legend.title=element_blank())  + scale_colour_viridis(discrete=TRUE) + xlab(xlab) + ylab(ylab)+  
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max) + facet_wrap(~factor(Interval))
  
  print(p9)
  dev.copy(png, paste(i, "8km_train_loocv_predperf_months_separately.png", sep="_"), width=900, height=750)
  dev.off()
  
  
  
  # Residual plot
  p10 <- ggplot(obspred, aes(x=obs, y=obs-pred)) +
    geom_point(shape = 16,  size=3) + geom_abline(slope = 0) + 
    
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    
    theme_pub  + theme(legend.title=element_blank())  + scale_colour_viridis(discrete=FALSE) + xlab("obs") + 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max) 
  
  print(p10)
  dev.copy(png, paste(i, "8km_train_loocv_predperf_residuals.png", sep="_"), width=900, height=750)
  dev.off()
  
  
  
  
  
  
  
  
  ### Visualize with bins instead
  # p11 <- ggplot(obspred) + geom_bin2d(aes(x=pred, y=obs),bins=20) + geom_abline(slope = 1) + 
  #   xlim(scale_min, scale_max) + ylim(scale_min, scale_max)  + theme_pub + 
  #   scale_fill_continuous(type = "viridis",  trans = 'reverse') + # could add trans='log' but then counts are harder to interpret?
  #   annotate(label = sprintf("\n MAE = %.1f \n Rsquared = %.2f \n RMSE = %.1f \n ", 
  #                            r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
  #                            r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
  #                            r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) 
  # print(p11)
  
  
  p11 <- ggplot(obspred) + geom_bin2d(aes(x=pred, y=obs),bins=50) + geom_abline(slope = 1) + 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)  + theme_pub + 
    scale_fill_viridis(direction=-1,  trans = 'log', labels = scales::number_format(accuracy = 1)) + # could add trans='log' but then counts are harder to interpret?
    annotate(label = sprintf("\n MAE = %.1f  \n MBE = %.1f    \n Rsquared = %.2f   \n RMSE = %.1f   \n ", 
                             r_stats  %>% dplyr::select(MAE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(MBE) %>% as.character() %>% as.numeric(), 
                             r_stats  %>% dplyr::select(Rsquared) %>% as.numeric(),
                             r_stats  %>% dplyr::select(RMSE) %>% as.numeric()), geom = "text", x = -Inf, y = Inf, size = 8, hjust = 0, vjust = 1) +
    labs(fill="Log counts") + xlab(xlab) + ylab(ylab)
  
  print(p11)
  dev.copy(png, paste(i, "8km_train_loocv_predperf_bins.png", sep="_"), width=500, height=400)
  dev.off()
  
  
  p12 <- ggplot(obspred) + geom_bin2d(aes(x=pred, y=obs),bins=50) + geom_abline(slope = 1) + 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)  + theme_pub + 
    scale_fill_viridis(direction=-1,  trans = 'log', labels = scales::number_format(accuracy = 1)) + # could add trans='log' but then counts are harder to interpret?
    labs(fill="Log counts") + facet_wrap(~Country) + xlab(xlab) + ylab(ylab)
  
  print(p12)
  dev.copy(png, paste(i, "8km_train_loocv_predperf_bins_country.png", sep="_"), width=900, height=750)
  dev.off()
  
  
  obspred$Disturbance2 <- as.character(obspred$Disturbance)
  obspred$Disturbance2 <- ifelse(obspred$Disturbance2=="Drainage", "NA/No", obspred$Disturbance2)
  p13 <- ggplot(obspred) + geom_bin2d(aes(x=pred, y=obs),bins=50) + geom_abline(slope = 1) + 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)  + theme_pub + 
    scale_fill_viridis(direction=-1,  trans = 'log', labels = scales::number_format(accuracy = 1)) + # could add trans='log' but then counts are harder to interpret?
    labs(fill="Log counts") + facet_wrap(~Disturbance2) + xlab(xlab) + ylab(ylab)
  
  print(p13)
  dev.copy(png, paste(i, "8km_train_loocv_predperf_bins_disturbance.png", sep="_"), width=900, height=750)
  dev.off()
  
  
  
  
  
  
  ### Model fit 
  library("randomForest")
  predtest <- predict(mod, obspred)
  obspred$predtest <- predtest   /1 * (maxval - minval) + minval
  p11 <- ggplot(obspred) + geom_bin2d(aes(x=predtest, y=obs), bins=30) + geom_abline(slope = 1) + 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max) + facet_wrap(~Country)
  
  print(p11)
  dev.copy(png, paste(i, "8km_train_loocv_modelfit_bins_country.png", sep="_"), width=900, height=750)
  dev.off()
  
  
  
  
  p11 <- ggplot(obspred) + geom_bin2d(aes(x=predtest, y=obs),bins=50) + geom_abline(slope = 1) + 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max)  + theme_pub + 
    scale_fill_viridis(direction=-1,  trans = 'log', labels = scales::number_format(accuracy = 1))+
    labs(fill="Log counts") + xlab(xlab) + ylab(ylab)
  
  print(p11)
  dev.copy(png, paste(i, "8km_train_loocv_modelfit_bins.png", sep="_"), width=500, height=400)
  dev.off()
  
  
  
  
  # if grouped to annual mean...
  obspredtest <- obspred %>% group_by(Study_ID_Short, Country, Biome, Meas_year) %>% summarize(obs_sum=sum(obs), predtest_sum=sum(predtest), n=n()) %>% filter(n==12)
  ggplot(obspredtest) + geom_bin2d(aes(x=predtest_sum, y=obs_sum), bins=30) + geom_abline(slope = 1) + 
    xlim(scale_min, scale_max) + ylim(scale_min, scale_max) + facet_wrap(~paste(Country, Biome))
  
  
  # Sweden and USA have quite many source observations that are predicted as larger Co2 sources
  obspredtest_long <- pivot_longer(obspredtest, cols=c(5,6))
  ggplot(obspredtest_long) + geom_density(aes(value, fill=name, color=name), alpha=0.5) + facet_wrap(~Country) + facet_wrap(~paste(Country, Biome)) + geom_vline(xintercept=0)
  
  
  
  
  
  ### Visualize with boxplots
  str(obspred)
  obspred_long <- pivot_longer(obspred, c("obs", "pred", "predtest"))
  obspred_long$name <- ifelse(obspred_long$name=="obs", "Observed", obspred_long$name)
  obspred_long$name <- ifelse(obspred_long$name=="pred", "Predicted (CV)", obspred_long$name)
  obspred_long$name <- ifelse(obspred_long$name=="predtest", "Predicted (no CV)", obspred_long$name)
  obspred_long$Biome_name <- paste(obspred_long$Biome, obspred_long$name)
  obspred_long$Biome_name <- factor(obspred_long$Biome_name, levels=c("Boreal Observed", "Boreal Predicted (no CV)", "Boreal Predicted (CV)",
                                                                      "Tundra Observed", "Tundra Predicted (no CV)", "Tundra Predicted (CV)"))
  
  
  
  # Define x and y lab titles for the plot
  if (i=="GPP_gC_m2") {
    ylab2 = expression(paste("GPP g C m"^{-2}, month^{-1}))
    
  }
  
  if (i=="NEE_gC_m2") {
    ylab2 = expression(paste("NEE g C m"^{-2}, month^{-1}))
    
  }
  
  if (i=="Reco_gC_m2") {
    ylab2 = expression(paste("Reco g C m"^{-2}, month^{-1}))
    
  }
  
  p14 <- ggplot(obspred_long) + geom_boxplot(aes(x=factor(Interval), y=value), col="black") + facet_wrap(~paste(Biome_name))  + theme_pub + ylab(ylab2) + xlab("Month")
  
  print(p14)
  dev.copy(png, paste(i, "8km_train_loocv_predperf_boxplot.png", sep="_"), width=900, height=750)
  dev.off()
  
  
  
  
  
 
  
  
  print("final pred perf figures done")
  
  
  
  ### Variable importance
  # #varimp <- caret::varImp(mod$finalModel, scale=TRUE)
  # all_varImp <- importance(mod$finalModel, scale=TRUE) %>% data.frame() ### TEMPORARY
  
  # extract varimp results - this should work
  all_varImp <- caret::varImp(mod$finalModel, scale=TRUE) %>% data.frame() ##### USE VIP HERE INSTEAD!!!!!!
  
  #https://stats.stackexchange.com/questions/109270/caret-varimp-for-randomforest-model
  
  
  # ## to use other functions of the package randomForest, convert class back
  # class(mod$finalModel) <- "randomForest"
  # importance(mod$finalModel, scale=TRUE) ## importance measure from the standard RF
  # 
  # 
  # class(mod) <- "randomForest"
  # importance(mod$finalModel, scale=TRUE)
  # 
  #   
  # mod$finalModel$importance[,1]/mod$finalModel$importanceSD
  # 
  # varImp(mod$finalModel)
  #   
  # library("randomForest")
  # importance(mod$finalModel, type=2, class=NULL, scale=TRUE)
  #   
  #   
  # varImpPlot(mod)### Variable importance calculation
  
  
  
  all_varImp$Variable2 <- row.names(all_varImp)
  
  all_varImp$Importance <- all_varImp$Overall
  
  
  # change the name
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="srad_terraclimate_sites", "Solar radiation", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="vpd_terraclimate_sites", "Vapor pressure deficit", all_varImp$Variable2)
  
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="tmean_terraclimate_sites", "Monthly air temperature", all_varImp$Variable2)
  
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="pr_terraclimate_sites", "Precipitation", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="pdsi_terraclimate_sites", "PDSI", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="tmean_TerraClimate_averages", "Mean annual air temperature over 1961-1990", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="ppt_TerraClimate_averages", "Mean annual precipitation over 1961-1990", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="trend_20yrprior_terra_change_id", "Rolling last 20-year air temperature trend", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="terra_trend_19601990", "Annual mean air temperature trend over 1961-1990", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="ndvi_trend_19812010", "June-August mean NDVI trend over 1982-2010", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Barrow_CO2_conc_Barrow_CO2conc", "Atmospheric CO2", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Snow.cover_era5_soilmoist_temp_snow", "Snow cover", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Snow.depth_era5_soilmoist_temp_snow", "Snow depth", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Soil.temperature.level.1_era5_soilmoist_temp_snow", "Soil temperature", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow", "Volumetric soil water", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="NDVI_whittaker_constant_monthly_mean", "NDVI", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality", "Land surface temperature (daytime)", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="water_ground_MCD43A4_annual_water_ground_sites_low_quality", "June-August average NDWI", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="water_vegetation_MCD43A4_annual_water_vegetation_sites_low_quality", "June-August average NDII", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m", "Compound topographic index", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="dtm_rough.scale_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m", "Topographic roughness scale", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="aboveground_biomass_carbon_2010_Above_belowground_biomass", "Aboveground biomass", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="belowground_biomass_carbon_2010_Above_belowground_biomass", "Belowground biomass", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Percent_NonTree_Vegetation_MOD44B_sites", "Percent non-tree vegetation", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Percent_NonVegetated_MOD44B_sites", "Percent non-vegetated", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Percent_Tree_Cover_MOD44B_sites", "Percent tree cover", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged", "Vegetation type", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="PHIHOX_M_sl1_250m_ll_SoilGrids", "Soil pH", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="BLDFIE_M_sl1_250m_ll_SoilGrids", "Soil bulk density", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="sol_watercontent.1500kPa_usda.3c2a1a_m_250m_b0..0cm_1950..2017_v0.1_SoilGrids_watercontent", "Soil water content", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="TKWP_Thermokarst", "Wetland thermokarst", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="TKHP_Thermokarst", "Hillslope thermokarst", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="forest_age_class_forest_age_sites", "Vegetation age", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH", "Permafrost probability", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="SoilGrids_SOC_SoilGrids_SOCstock", "Soil organic carbon stock", all_varImp$Variable2)
  
  
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Percent_TreeCover_AVHRR_VCF5KYR", "Percent tree cover", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Percent_NonTree_Vegetation_AVHRR_VCF5KYR", "Percent non-tree vegetation", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="Percent_NonVegetated_AVHRR_VCF5KYR", "Percent non-vegetated", all_varImp$Variable2)
  all_varImp$Variable2 <- ifelse(all_varImp$Variable2=="ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled", "NDVI", all_varImp$Variable2)
  
  all_varImp <- all_varImp[order(all_varImp$Importance, decreasing=TRUE),]
  all_varImp$Variable2 <- factor(all_varImp$Variable2, levels=unique(all_varImp$Variable2))
  
  
  p1 <- ggplot(all_varImp) + geom_bar(aes(x=Variable2, y=Importance), stat="identity")   + 
    coord_flip() + scale_x_discrete(limits = rev(levels(all_varImp$Variable2)))+ 
    theme_pub #+ ggtitle(title2)
  
  # Print out
  setwd("/home/master/abcfluxv1_modeling/figures/")
  print(p1)
  
  dev.copy(png, paste( i, "qrf_vip.png", sep="_"), width=1300, height=1100)
  dev.off()
  
  
  # EPS
  
  # Print out
  setwd("/home/master/abcfluxv1_modeling/figures/")
  ggsave(paste( i,  "8km_qrf_vip.eps", sep="_"), device=cairo_ps, p1, width=35, height=23, units=c("cm"))
  
  
  write.csv(all_varImp, paste("/home/master/abcflux_modeling/results/", i, "8km_qrf_vip.csv", sep="_"), row.names=FALSE)
  
  
      
     
  
}




####### TEMPORARILY UNCOMMENTED

### Partial dependence plots need to be done in a separate loop



for (i in resp_vars) {


    # i <- "NEE_gC_m2"
    # i <- "GPP_gC_m2"
    # i <- "Reco_gC_m2"


    # Add data
    
    modeldata2 <- d[,c("Study_ID_Short", "id", i, Baseline_vars_20km)]
    modeldata2 <- na.omit(modeldata2) 
    
    # create a row ID
    modeldata2$samplerow <- seq(1, length(modeldata2$Study_ID_Short), by=1)
    
    
    # merge back other information
    modeldata22 <- merge(modeldata2[ , !(names(modeldata2) %in%  c("Study_ID_Short",  i, Baseline_vars_20km))], d, by="id")
    
    
    
    # i <- "NEE_gC_m2"
    # i <- "GPP_gC_m2"
    # i <- "Reco_gC_m2" 
    
    
    
    
    print(paste0("/home/master/abcfluxv1_modeling/results/", paste(i,  "8km_qrf", "train_loocv_norm",sep="_"), ".rds"))
    # Load model files
    mod <- readRDS(paste0("/home/master/abcfluxv1_modeling/results/", paste(i, "8km_qrf",  "train_loocv_norm",  sep="_"), ".rds"))
    mod2 = mod$finalModel #pull out the quantile regression object
    
    
    ### Predictive performance plots ###
    
    
    # Define x and y lab titles for the plot
    
    
    if (i=="GPP_gC_m2") {
      ylab = expression(paste("Observed GPP g C m"^{-2}, month^{-1}))
      xlab = expression(paste("Predicted GPP g C m"^{-2}, month^{-1}, " (8 km model)"))
      
    }
    
    if (i=="NEE_gC_m2") {
      ylab = expression(paste("Observed NEE g C m"^{-2}, month^{-1}))
      xlab = expression(paste("Predicted NEE g C m"^{-2}, month^{-1}, " (8 km model)"))
      
    }
    
    if (i=="Reco_gC_m2") {
      ylab = expression(paste("Observed Reco g C m"^{-2}, month^{-1}))
      xlab = expression(paste("Predicted Reco g C m"^{-2}, month^{-1}, " (8 km model)"))
      
    }
    
    
    
    # # Extract the final model details
    preds <- mod$pred %>%
      data.frame()
    
    
    
    
    # Merge
    
    obspred <- merge(modeldata22, preds, by.x="samplerow", by.y="rowIndex")
    # plot(obspred$NEE_gC_m2_scaled, obspred$obs) # yes, identical
    # cor.test(obspred$NEE_gC_m2, obspred$obs)
    
    
    nvars <- length(Baseline_vars_20km)


    library("randomForest")

    for (nvar in 1:nvars) {

      pd1 <- partial(mod, pred.var = Baseline_vars_20km[nvar], train=obspred)  # don't set plot = TRUE

      rug_data <- obspred[Baseline_vars_20km[nvar]]

      # variable names
      # change the name
      # change the name
      Baseline_vars_20km_name <- Baseline_vars_20km
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="srad_terraclimate_sites", "Solar radiation/10 W m-2	", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="vpd_terraclimate_sites", "Vapor pressure deficit/100 kPa", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="pr_terraclimate_sites", "Precipitation mm", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="tmean_terraclimate_sites", "Air temperature °C", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="pdsi_terraclimate_sites", "PDSI/100", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="tmean_TerraClimate_averages", "Mean annual air temperature over 1961-1990 ?C", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="ppt_TerraClimate_averages", "Mean annual precipitation over 1961-1990 mm", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="trend_20yrprior_terra_change_id", "Rolling last 20-year air temperature trend/10 ?C", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="terra_trend_19601990", "Annual mean air temperature trend over 1961-1990/10 ?C", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="ndvi_trend_19812010", "June-August mean NDVI trend over 1982-2010", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Barrow_CO2_conc_Barrow_CO2conc", "Atmospheric CO2", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Snow.cover_era5_soilmoist_temp_snow", "Snow cover %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Snow.depth_era5_soilmoist_temp_snow", "Snow depth m", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Soil.temperature.level.1_era5_soilmoist_temp_snow", "Soil temperature Kelvin", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Volumetric.soil.water.layer.1_era5_soilmoist_temp_snow", "Volumetric soil water m3 m-3", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="NDVI_whittaker_constant_monthly_mean", "NDVI", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="LST_Day_1km_MOD11A2v006_LST_Day_sites_low_quality", "Land surface temperature (daytime) Kelvin", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="water_ground_MCD43A4_annual_water_ground_sites_low_quality", "June-August average NDWI", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="water_vegetation_MCD43A4_annual_water_vegetation_sites_low_quality", "June-August average NDII", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="dtm_cti_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m", "Compound topographic index", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="dtm_rough.scale_merit.dem_m_250m_s0..0cm_2018_v1.0_MERIT_topo_indices_250m", "Topographic roughness scale", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="aboveground_biomass_carbon_2010_Above_belowground_biomass", "Aboveground biomass MgC ha-1", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="belowground_biomass_carbon_2010_Above_belowground_biomass", "Belowground biomass MgC ha-1", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Percent_NonTree_Vegetation_MOD44B_sites", "Percent non-tree vegetation %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Percent_NonVegetated_MOD44B_sites", "Percent non-vegetated %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Percent_Tree_Cover_MOD44B_sites", "Percent tree cover %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged", "Vegetation type", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="PHIHOX_M_sl1_250m_ll_SoilGrids", "Soil pH", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="BLDFIE_M_sl1_250m_ll_SoilGrids", "Soil bulk density kg m-3", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="sol_watercontent.1500kPa_usda.3c2a1a_m_250m_b0..0cm_1950..2017_v0.1_SoilGrids_watercontent", "Soil water content %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Number_of_days_since_fire_classes_MCD64A1_sites_cleaned", "Time since fire (no burn - old burn)", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="TKWP_Thermokarst", "Wetland thermokarst (low - high coverage)", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="TKHP_Thermokarst", "Hillslope thermokarst (low - high coverage)", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="forest_age_class_forest_age_sites", "Vegetation age (old - young)", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH_UiO_PEX_20181128_2000_2016_NH", "Permafrost probability %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="SoilGrids_SOC_SoilGrids_SOCstock", "Soil organic carbon stock Tonnes ha-1", Baseline_vars_20km_name)
      
      
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Percent_TreeCover_AVHRR_VCF5KYR", "Percent tree cover %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Percent_NonTree_Vegetation_AVHRR_VCF5KYR", "Percent non-tree vegetation %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="Percent_NonVegetated_AVHRR_VCF5KYR", "Percent non-vegetated %", Baseline_vars_20km_name)
      Baseline_vars_20km_name <- ifelse(Baseline_vars_20km_name=="ndvi3g_lowest_gapfilled_mean_GIMMS3g_NDVI_sites_low_quality_gapfilled", "NDVI", Baseline_vars_20km_name)

      if (!is.factor(obspred[, Baseline_vars_20km[nvar]])) {

        pdp_plot1 <- ggplot(pd1, aes(x=pd1[, 1], y=pd1[, 2])) +
          geom_line(size=1) + theme_pub  +labs(y="yhat", x=Baseline_vars_20km_name[nvar]) +
          theme(legend.position = "none") +
          geom_rug(data=rug_data, aes(x=rug_data[, 1]), inherit.aes = FALSE) #+ geom_rug(col=rgb(.5,0,0,alpha=.2)) # + ggtitle(title)
      } else if  (is.factor(obspred[, Baseline_vars_20km[nvar]])) {

        pdp_plot1 <- ggplot(pd1, aes(x=pd1[, 1], y=pd1[, 2])) +
          geom_point(size=3) + theme_pub  +labs(y="yhat", x=Baseline_vars_20km_name[nvar]) +
          theme(legend.position = "none") # + ggtitle(title)
      }


      dev.copy(png, paste( i, Baseline_vars_20km[nvar],  "8km_qrf_pdp.png", sep="_"), width=600, height=400)
      print(pdp_plot1)
      dev.off()




    } # nvars loop




} # resp vars loop
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # # Partial dependence plots - testing
# # 
# # # First baseline var
# # pd1 <- partial(gbmFit, pred.var = Baseline_vars_1km[36])  # don't set plot = TRUE
# # pd2 <- partial(rfFit, pred.var = Baseline_vars_1km[36])
# # pd3 <- partial(svmFit, pred.var = Baseline_vars_1km[36]) # , train=d
# # 
# # # ggplot2
# # pd1$Model <- "GBM"  # add new column
# # pd2$Model <- "RF"
# # pd3$Model <- "SVM"
# # 
# # pd.all1 <- rbind(pd1, pd2, pd3)  # bind rows
# # 
# # 
# # if(km=="1km") {
# #   
# #   # First baseline var
# #   pd1 <- partial(gbmFit, pred.var = Baseline_vars_1km[nvar])  # don't set plot = TRUE
# #   pd2 <- partial(rfFit, pred.var = Baseline_vars_1km[nvar])
# #   pd3 <- partial(svmFit, pred.var = Baseline_vars_1km[nvar]) # , train=d
# #   
# #   # ggplot2
# #   pd1$Model <- "GBM"  # add new column
# #   pd2$Model <- "RF"
# #   pd3$Model <- "SVM"
# #   
# #   pdp_plot1 <- ggplot(pd.all1, aes(x=pd.all1[, 1], y=pd.all1[, 2], color = Model)) +
# #     geom_line(size=1) + theme_pub  +labs(y="yhat", x=Baseline_vars_1km[1]) + 
# #     scale_color_manual(values = c(viridis(3)[3],viridis(3)[2], viridis(3)[1]), guide = guide_legend(reverse=TRUE)) +
# #     theme(legend.position = "none") # + ggtitle(title)
# # } if (km=="20km") {
# #   
# #   # First baseline var
# #   pd1 <- partial(gbmFit, pred.var = Baseline_vars_20km[nvar])  # don't set plot = TRUE
# #   pd2 <- partial(rfFit, pred.var = Baseline_vars_20km[nvar])
# #   pd3 <- partial(svmFit, pred.var = Baseline_vars_20km[nvar]) # , train=d
# #   
# #   # ggplot2
# #   pd1$Model <- "GBM"  # add new column
# #   pd2$Model <- "RF"
# #   pd3$Model <- "SVM"
# #   
# #   pdp_plot1 <- ggplot(pd.all1, aes(x=pd.all1[, 1], y=pd.all1[, 2], color = Model)) +
# #     geom_line(size=1) + theme_pub  +labs(y="yhat", x=Baseline_vars_20km[1]) + 
# #     scale_color_manual(values = c(viridis(3)[3],viridis(3)[2], viridis(3)[1]), guide = guide_legend(reverse=TRUE)) +
# #     theme(legend.position = "none") # + ggtitle(title)
# #   
# # }
# # 
# # 
# # print(pdp_plot1)
# # 
# # dev.copy(png, paste0( i, km, nvar, "pdp.png", sep="_"), width=650, height=400)
# # dev.off()
# 
# 
# 
