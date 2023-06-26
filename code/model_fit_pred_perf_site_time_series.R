





.libPaths("/home/master/R/x86_64-pc-linux-gnu-library/4.2")

library("viridis")
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
d <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg_allsites_1km.csv")
d <- d %>% mutate_if(is.integer, as.numeric)

# Reverse gpp first so that GPP and Reco have both positive values
d$GPP_gC_m2 <- -d$GPP_gC_m2 

# Simplifying the land cover classes (for simplicity and because some classes have very limited amount of data)
d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==1, 31, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)

d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==90, 30, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)

d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- ifelse(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged==41, 160, d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)







# merge some site-level data back (only done in the loop - this can change the data order)
d_orig <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg.csv")
d_orig <- subset(d_orig, select=c("Study_ID_Short", "Biome", "Disturbance", "Country")) %>% unique()
duplicated(d_orig$Study_ID_Short)

# test <- merge(d, d_orig, by="Study_ID_Short")
# unique(d$Study_ID_Short) %in% unique(test$Study_ID_Short)
# unique(test$Study_ID_Short)




d$Study_ID_Short[d$Study_ID_Short=="L\xf3pez-Blanco_GL-NuF_tower1"] <- "Lopez-Blanco_GL-NuF_tower1"
d$Study_ID_Short[d$Study_ID_Short=="L\xf3pez-Blanco_GL-ZaF_tower1"] <- "Lopez-Blanco_GL-Zaf_tower1"

# Define the factors
d$Study_ID_Short <- factor(d$Study_ID_Short)
d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged <- as.factor(d$ESACCI_cavm_general_ESAwaterfix_broadevfix_mixfix_cropfix_nowaterglacier_ESACCI_CAVM_merged)
d$Interval <- as.factor(d$Interval)





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
resp_vars <- c("NEE_gC_m2",  "GPP_gC_m2", "Reco_gC_m2") # 





### Set folder for results
setwd("/home/master/abcfluxv1_modeling/figures")



# For figures

theme_pub <- theme_bw() + theme(panel.border=element_rect(size=1, colour="black"),
                                axis.text=element_text(size=20, face="bold"),
                                plot.title = element_text(size = 20, face = "bold"),
                                axis.title=element_text(size=20, face="bold"), 
                                plot.subtitle=element_text(size=20, face="bold", color="black"), 
                                strip.text.x = element_text(size = 20, face="bold"),
                                legend.text=element_text(size=20, face="bold"), legend.title=element_text(size=20))




theme_pub <- theme_bw() + theme(panel.border=element_rect(size=1, colour="black"),
                                axis.text=element_text(size=20, face="bold"),
                                plot.title = element_text(size = 20, face = "bold"),
                                axis.title=element_text(size=20, face="bold"), 
                                plot.subtitle=element_text(size=20, face="bold", color="black"), 
                                strip.text.x = element_text(size = 20, face="bold"),
                                legend.text=element_text(size=20, face="bold"), legend.title=element_text(size=20))

# function to get the flux model data
flux_model <- function(i) {
  
  #i <- "NEE_gC_m2"
  #i <- "GPP_gC_m2"
  #i <- "Reco_gC_m2"
  
  # Add data
  
  modeldata2 <- d[,c("Study_ID_Short", "id", i, Baseline_vars_1km)]
  modeldata2 <- na.omit(modeldata2) 
  
  # create a row ID
  modeldata2$samplerow <- seq(1, length(modeldata2$Study_ID_Short), by=1)
  
  
  # merge back other information
  modeldata22 <- merge(modeldata2[ , !(names(modeldata2) %in%  c("Study_ID_Short",  i, Baseline_vars_1km))], d, by="id")
  
  
  
  # i <- "NEE_gC_m2"
  # i <- "GPP_gC_m2"
  # i <- "Reco_gC_m2" 
  
  
  
  
  print(paste0("/home/master/abcfluxv1_modeling/results/", paste(i,  "1km_rf", "train_loocv",sep="_"), ".rds"))
  # Load model files
  mod <- readRDS(paste0("/home/master/abcfluxv1_modeling/results/", paste(i, "1km_rf",  "train_loocv",  sep="_"), ".rds"))
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
  
  
  # Model fit
  modeldata22$pred_fit <- predict(mod2, modeldata22[, Baseline_vars_1km])
  
  # Merge
  obspred <- merge(modeldata22, preds, by.x="samplerow", by.y="rowIndex")
  plot(obspred[, i], obspred$obs) # yes, identical - nope, not always!
  # cor.test(obspred$NEE_gC_m2, obspred$obs)
  # 
  
  # 
  # obspred <- merge(obspred, d_orig, by="Study_ID_Short")
  
  # dates
  obspred$Date <- as.Date(paste(obspred$Meas_year, obspred$Interval, "01", sep="-"))
  
  # wrangling
  obspred_long <- obspred %>% pivot_longer(cols=c("obs", "pred", "pred_fit"), names_to="prediction_type", values_to="predictions")

  obspred_long$prediction_type <- ifelse(obspred_long$prediction_type=="obs", "Observation", obspred_long$prediction_type)
  obspred_long$prediction_type <- ifelse(obspred_long$prediction_type=="pred", "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)", obspred_long$prediction_type)
  obspred_long$prediction_type <- ifelse(obspred_long$prediction_type=="pred_fit", "Prediction based on \nthe full model training data \n(i.e., model fit)", obspred_long$prediction_type)
  
  obspred_long$prediction_type <- factor(obspred_long$prediction_type,   levels = c("Observation", "Prediction based on \nthe full model training data \n(i.e., model fit)", "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)")
)
  #levels(obspred_long$prediction_type) <- c("Observation", "Prediction based on \nthe full model training data \n(i.e., model fit)", "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)")
  
}


### NEE
flux_model("NEE_gC_m2")
i <- "NEE_gC_m2"

### EML time series. Observations in black, model fit in red, pred perf in grey

obspred_long_sub <- subset(obspred_long, obspred_long$Study_ID_Short=="Schuur_US-EML_tower1")

ggplot(data=obspred_long_sub) + geom_line(aes(x=Date, y=predictions, col=prediction_type), size=1.5, alpha=0.5) + 
  theme_pub + scale_color_viridis(discrete = TRUE, option = "D") + ylab(expression(bold(paste("NEE g C m"^{-2}, month^{-1})))) + 
  ggtitle("US-EML") +labs(colour="")

setwd("/home/master/abcfluxv1_modeling/code")
dev.copy(png, paste0("../figures/", i, "_model_fit_pred_perf_eml.png"), width=850, height=400)
dev.off()

### Cherski

obspred_long_sub <- subset(obspred_long, obspred_long$Study_ID_Short=="Goeckede_RU-Ch2_tower2") 
# calculate all days
all_days = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days$prediction_type <- "Observation"
all_days2 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days2$prediction_type <- "Prediction based on \nthe full model training data \n(i.e., model fit)"
all_days3 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days3$prediction_type <- "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)"
all_days <- rbind(all_days, all_days2, all_days3)
# join to original data
library(dplyr)
obspred_long_sub = left_join(all_days, obspred_long_sub, by = c("Date", "prediction_type"))
obspred_long_sub$prediction_type <- factor(obspred_long_sub$prediction_type,   levels = c("Observation", "Prediction based on \nthe full model training data \n(i.e., model fit)", "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)")
)
ggplot(data=obspred_long_sub) + geom_line(aes(x=Date, y=predictions, col=prediction_type), size=1.5, alpha=0.5) + 
  theme_pub + scale_color_viridis(discrete = TRUE, option = "D") + ylab(expression(bold(paste("NEE g C m"^{-2}, month^{-1})))) + 
  ggtitle("RU-Ch2") +labs(colour="")
dev.copy(png, paste0("../figures/", i, "_model_fit_pred_perf_ruch2.png"), width=850, height=400)
dev.off()

### Hyytiälä

obspred_long_sub <- subset(obspred_long, obspred_long$Study_ID_Short=="Vesala_FI-Hyy_tower1") 

ggplot(data=obspred_long_sub) + geom_line(aes(x=Date, y=predictions, col=prediction_type), size=1.5, alpha=0.5) + 
  theme_pub + scale_color_viridis(discrete = TRUE, option = "D") + ylab(expression(bold(paste("NEE g C m"^{-2}, month^{-1})))) + 
  ggtitle("FI-Hyy") +labs(colour="")

dev.copy(png, paste0("../figures/", i, "_model_fit_pred_perf_fihyy.png"), width=850, height=400)
dev.off()


### Lopez-Blanco_GL-Zaf_tower1

obspred_long_sub <- subset(obspred_long, obspred_long$Study_ID_Short=="Lopez-Blanco_GL-Zaf_tower1") 
# calculate all days
all_days = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days$prediction_type <- "Observation"
all_days2 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days2$prediction_type <- "Prediction based on \nthe full model training data \n(i.e., model fit)"
all_days3 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days3$prediction_type <- "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)"
all_days <- rbind(all_days, all_days2, all_days3)
# join to original data
library(dplyr)
obspred_long_sub = left_join(all_days, obspred_long_sub, by = c("Date", "prediction_type"))
obspred_long_sub$prediction_type <- factor(obspred_long_sub$prediction_type,   levels = c("Observation", "Prediction based on \nthe full model training data \n(i.e., model fit)", "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)")
)
ggplot(data=obspred_long_sub) + geom_line(aes(x=Date, y=predictions, col=prediction_type), size=1.5, alpha=0.5) + 
  theme_pub + scale_color_viridis(discrete = TRUE, option = "D") + ylab(expression(bold(paste("NEE g C m"^{-2}, month^{-1})))) + 
  ggtitle("GL-Zaf") +labs(colour="")
dev.copy(png, paste0("../figures/", i, "_model_fit_pred_perf_glzaf.png"), width=850, height=400)
dev.off()





# Maximov_RU-SkP_tower1
obspred_long_sub <- subset(obspred_long, obspred_long$Study_ID_Short=="Maximov_RU-SkP_tower1" & Meas_year> 2005) 
# calculate all days
all_days = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days$prediction_type <- "Observation"
all_days2 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days2$prediction_type <- "Prediction based on \nthe full model training data \n(i.e., model fit)"
all_days3 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days3$prediction_type <- "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)"
all_days <- rbind(all_days, all_days2, all_days3)
# join to original data
library(dplyr)
obspred_long_sub = left_join(all_days, obspred_long_sub, by = c("Date", "prediction_type"))
obspred_long_sub$prediction_type <- factor(obspred_long_sub$prediction_type,   levels = c("Observation", "Prediction based on \nthe full model training data \n(i.e., model fit)", "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)")
)
ggplot(data=obspred_long_sub) + geom_line(aes(x=Date, y=predictions, col=prediction_type), size=1.5, alpha=0.5) + 
  theme_pub + scale_color_viridis(discrete = TRUE, option = "D") + ylab(expression(bold(paste("NEE g C m"^{-2}, month^{-1})))) + 
  ggtitle("RU-SkP") +labs(colour="")
dev.copy(png, paste0("../figures/", i, "_model_fit_pred_perf_ruskp.png"), width=850, height=400)
dev.off()



#Helbig_CA-SCB_tower1
obspred_long_sub <- subset(obspred_long, obspred_long$Study_ID_Short=="Helbig_CA-SCB_tower1" & Meas_year> 2005) 
# calculate all days
all_days = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days$prediction_type <- "Observation"
all_days2 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days2$prediction_type <- "Prediction based on \nthe full model training data \n(i.e., model fit)"
all_days3 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
all_days3$prediction_type <- "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)"
all_days <- rbind(all_days, all_days2, all_days3)
# join to original data
library(dplyr)
obspred_long_sub = left_join(all_days, obspred_long_sub, by = c("Date", "prediction_type"))
obspred_long_sub$prediction_type <- factor(obspred_long_sub$prediction_type,   levels = c("Observation", "Prediction based on \nthe full model training data \n(i.e., model fit)", "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)")
)
ggplot(data=obspred_long_sub) + geom_line(aes(x=Date, y=predictions, col=prediction_type), size=1.5, alpha=0.5) + 
  theme_pub + scale_color_viridis(discrete = TRUE, option = "D") + ylab(expression(bold(paste("NEE g C m"^{-2}, month^{-1})))) + 
  ggtitle("CA-SCB") +labs(colour="")
dev.copy(png, paste0("../figures/", i, "_model_fit_pred_perf_cascb.png"), width=850, height=400)
dev.off()




# # Maximov_RU-Elg_tower1
# 
# obspred_long_sub <- subset(obspred_long, obspred_long$Study_ID_Short=="Maximov_RU-Elg_tower1" & Meas_year> 2005) 
# # calculate all days
# all_days = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
# all_days$prediction_type <- "Observation"
# all_days2 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
# all_days2$prediction_type <- "Prediction based on \nthe full model training data \n(i.e., model fit)"
# all_days3 = data.frame(Date = seq.Date(from = min(obspred_long_sub$Date), to = max(obspred_long_sub$Date), by = "month"))
# all_days3$prediction_type <- "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)"
# all_days <- rbind(all_days, all_days2, all_days3)
# # join to original data
# library(dplyr)
# obspred_long_sub = left_join(all_days, obspred_long_sub, by = c("Date", "prediction_type"))
# obspred_long_sub$prediction_type <- factor(obspred_long_sub$prediction_type,   levels = c("Observation", "Prediction based on \nthe full model training data \n(i.e., model fit)", "Prediction based on \nmodel training data \nexcluding the site \n(i.e., model predictive perf.\n based on CV)")
# )
# ggplot(data=obspred_long_sub) + geom_line(aes(x=Date, y=predictions, col=prediction_type), size=1.5, alpha=0.5) + 
#   theme_pub + scale_color_viridis(discrete = TRUE, option = "D") + ylab(expression(bold(paste("NEE g C m"^{-2}, month^{-1})))) + 
#   ggtitle("RU-SkP") +labs(colour="")
