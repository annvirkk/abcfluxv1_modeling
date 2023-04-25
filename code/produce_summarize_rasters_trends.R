
.libPaths("/home/master/R/x86_64-pc-linux-gnu-library/4.2")

library("stringr")
library("stringi")
library("terra")
library(dplyr)
library(googleCloudStorageR)
library(purrr)
library("zyp")

setwd("/home/master/")
terraOptions(memfrac=0.9, tempdir = "temp/") 

# Google cloud settings
my_project_id <- "top-operand-328213"
gcs_list_buckets(my_project_id)
gcs_global_bucket("abcflux_modeling_files")
contents <- gcs_list_objects()
gcs_upload_set_limit(50000000L) # increasing data size limit for transferring data to google cloud

### 8 km predictions ###

## In situ data
d <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg_allsites.csv")
d <- d %>% mutate_if(is.integer, as.numeric)
# reverse gpp first, for simplicity
d$GPP_gC_m2 <- -d$GPP_gC_m2 


# regions

eco <- vect("/home/master/masking_summary_rasters/Ecoregions2017_tundraboreal.shp")
eco <- aggregate(eco, "BIOME_NUM")
crs(eco) <- "epsg:4326"
eco <- project(eco, "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" )

region <- vect("/home/master/masking_summary_rasters/combined_regions_TreatVirkkala_v2.shp")
region <- project(region, "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" )









### calculate average flux rasters ###

avg_rasters <- function(flux, pattern1, folder) {



  setwd("/home/master/")
  files_to_download <- grep("*.csv", contents$name, value = TRUE)
  files_to_download2 <- files_to_download[grepl(flux, files_to_download)]
  files_to_download2 <- files_to_download2[grepl(paste0("predictions_8km/csv/", folder, "/"), files_to_download2)]
  files_to_download2 <- files_to_download2[!(grepl("train_loocv", files_to_download2)) ]


  map(files_to_download2, function(x) gcs_get_object(x, saveToDisk = x, overwrite = TRUE))

  
  
  
  
  # YEARLY LOOP STARTS HERE
  
  for (y in pattern1) {
    
    
    
    print(y)
    ## temporal splits

    # read to a data frame
    setwd(paste0("/home/master/predictions_8km/csv/", folder, "/"))
    files <- list.files(, pattern=flux)
    files <- files[!(grepl("train_loocv", files)) ]
    df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
    
    
    # rename
    columns <- NA
    for (f in files) {
      column <- c("x", "y", f)
      columns <- c(columns, column)
      
    }
    
    columns <- columns[2:length(columns)]
    names(df) <- columns
    
    
    # normalization values back
    #/ 10000 * (maxval - minval) + minval
    # note that the predictions have still been multiplied by 1000
    
    
    minval <- ifelse(flux=="NEE_gC_m2", min(d$NEE_gC_m2, na.rm=TRUE), NA)
    minval <- ifelse(flux=="GPP_gC_m2", min(d$GPP_gC_m2, na.rm=TRUE), minval)
    minval <- ifelse(flux=="Reco_gC_m2", min(d$Reco_gC_m2, na.rm=TRUE), minval)
    
    
    maxval <- ifelse(flux=="NEE_gC_m2", max(d$NEE_gC_m2, na.rm=TRUE), NA)
    maxval <- ifelse(flux=="GPP_gC_m2", max(d$GPP_gC_m2, na.rm=TRUE), maxval)
    maxval <- ifelse(flux=="Reco_gC_m2", max(d$Reco_gC_m2, na.rm=TRUE), maxval)
    
    # normalization values back
    # every third column, do this
    cols <- seq(3, length(files)*3, by=3)
    df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
    
    
    
    
    
    # months #
    months <- c("_01", "_02", "_03", "_04", "_05", "_06",
                "_07", "_08", "_09", "_10", "_11", "_12")
    
    for (m in months) {
      # m <- "_01"
      selected_columns <- str_detect(columns, m)
      selected_columns2 <- str_detect(columns, y)
      selected_columns<- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)

      selected_columns[1:2] <- TRUE # add coords
      df2 <- df[, selected_columns]
      print(colnames(df2)[3:length(colnames(df2))])
      
      df3 <- df2  # %>% mutate(mean = rowMeans(.[, 3:ncol(df2)])) # this would do annual averages
      df4 <- df3[, c(1:2, ncol(df3))]
      
      df4[, 3] <- unname(df4[, 3])
      
      rasterdf  <- as.matrix(df4[])
      r <- rast(rasterdf[], type="xyz")
      # plot(r)
      crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
      writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", m, "_", y, ".tif"), overwrite=TRUE, datatype="INT4S")
      
    }
    
    print("months done")
    
    # seasons: spring #
    months <- c("_03", "_04", "_05")
    
    # y <- "1982"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
    
    
    # y <- "1982"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
    
    
    
    # seasons: summer #
    months <- c("_06", "_07", "_08")
    
    
    # y <- "1982"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
    
    
    
    # y <- "1982"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
    
    
    
    
    
    # seasons: autumn #
    months <- c("_09", "_10", "_11")
    
    
    # y <- "2000"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
    
    
    # y <- "2000"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
    
    
    
    
    
    print("most of 3-month seasons done")
    
    
    
    # Calculate longer seasons (GS)
    
    # seasons: autumn #
    months <- c("_05","_06", "_07", "_08", "_09")
    
    
    # y <- "1982"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
    
    
    # y <- "1982"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
    
    print("GS done")
    
    
    
    
    # annual
    months <- c("_01", "_02", "_03", "_04", "_05", "_06",
                "_07", "_08", "_09", "_10", "_11", "_12")
    
    
    
    # y <- "1982"
    selected_columns <- str_detect(columns, paste(months, collapse = "|"))
    
    selected_columns2 <- str_detect(columns, y)
    
    selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
    
    selected_columns[1:2] <- TRUE # add coords
    df2 <- df[, selected_columns]
    print(colnames(df2)[3:length(colnames(df2))])
    
    df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
    
    df4 <- df3[, c(1, 2, ncol(df3))]
    
    rasterdf  <- as.matrix(df4[])
    r <- rast(rasterdf[], type="xyz")
    # plot(r)
    crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
    writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", "_", y, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
    
    
    
    print("annual done")
    
    
    
    # seasons: winter #
    # trickier because need to take previous year into account!
    months <- c("_01", "_02", "_12")
    
    
    
    # list and select files
    # files <- list.files(paste0("/home/master/predictions_8km/csv/", folder, "/"), pattern=paste0(pattern1, collapse="|")) # old way?
    files <- list.files(paste0("/home/master/predictions_8km/csv/", folder, "/"), pattern=y)
    files <- files[str_detect(files, flux)]
    files <- files[!str_detect(files, "_12")] # remove december from the same year
    files2 <- list.files(paste0("/home/master/predictions_8km/csv/", folder, "/"), pattern=paste0(as.numeric(y)-1, collapse="|"))
    files2 <- files2[str_detect(files2, "_12")]
    files2 <- files2[str_detect(files2, flux)]
    files <- c(files, files2)
    files <- files[str_detect(files, flux)]
    setwd(paste0("/home/master/predictions_8km/csv/", folder, "/"))
    
    # read to a data frame
    df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
    
    
    # rename
    columns <- NA
    for (f in files) {
      column <- c("x", "y", f)
      columns <- c(columns, column)
      
    }
    
    columns <- columns[2:length(columns)]
    names(df) <- columns
    
    
    
    # normalization values back
    # every third column, do this
    cols <- seq(3, length(files)*3, by=3)
    df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
    
    # skip 1982
    if (y!="1990") {
      
      # y <- "1983"
      # months 1-2
      selected_columns <- str_detect(columns, paste(months, collapse = "|"))

      selected_columns[1:2] <- TRUE # add coords
      df2 <- df[, selected_columns]
      print(colnames(df2)[3:length(colnames(df2))])
      
      df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
      
      df4 <- df3[, c(1, 2, ncol(df3))]
      
      rasterdf  <- as.matrix(df4[])
      r <- rast(rasterdf[], type="xyz")
      # plot(r)
      crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
      writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
      
      
      
      df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
      
      df4 <- df3[, c(1, 2, ncol(df3))]
      
      rasterdf  <- as.matrix(df4[])
      r <- rast(rasterdf[], type="xyz")
      # plot(r)
      crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
      writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
      
    }
    
    
    
    
    # seasons: non-growing season #
    # trickier because need to take previous year into account!
    months <- c("_01", "_02", "_03", "_04","_10","_11","_12")
    
    
    # list and select files
    files <- list.files(paste0("/home/master/predictions_8km/csv/", folder, "/"), pattern=paste0(y, collapse="|"))
    files <- files[!(str_detect(files, "_12") | str_detect(files, "_11") | str_detect(files, "_10"))]# remove december from the same year
    files2 <- list.files(paste0("/home/master/predictions_8km/csv/", folder, "/"), pattern=paste0(as.numeric(y)-1, collapse="|"))
    files2 <- files2[str_detect(files2, "_12") | str_detect(files2, "_11") | str_detect(files2, "_10") ]
    files <- c(files, files2)
    files <- files[str_detect(files, flux)]
    setwd(paste0("/home/master/predictions_8km/csv/", folder, "/"))
    
    # read to a data frame
    df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
    
    # rename
    columns <- NA
    for (f in files) {
      column <- c("x", "y", f)
      columns <- c(columns, column)
      
    }
    
    columns <- columns[2:length(columns)]
    names(df) <- columns
    
    
    
    
    # normalization values back
    # every third column, do this
    cols <- seq(3, length(files)*3, by=3)
    df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
    

      
      # skip 1982
      if (y!=1990) {
        
        # y <- "1983"
        # months 1-2
        selected_columns <- str_detect(columns, paste(months, collapse = "|"))
        
        #selected_columns2 <- str_detect(columns, y)
        
        
        #selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE | selected_columnsb==TRUE & selected_columnsb2==TRUE, TRUE, FALSE)
        
        selected_columns[1:2] <- TRUE # add coords
        df2 <- df[, selected_columns]
        print(colnames(df2)[3:length(colnames(df2))])
        
        df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
        
        df4 <- df3[, c(1, 2, ncol(df3))]
        
        rasterdf  <- as.matrix(df4[])
        r <- rast(rasterdf[], type="xyz")
        # plot(r)
        crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
        writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
        
        
        df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
        
        df4 <- df3[, c(1, 2, ncol(df3))]
        
        rasterdf  <- as.matrix(df4[])
        r <- rast(rasterdf[], type="xyz")
        # plot(r)
        crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
        writeRaster(r*100, paste0("/home/master/predictions_8km/raster/", folder, "/", flux, "_8km", stri_paste(months, collapse=''), "_", y, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
        
      }
      

    
    
    

    
    
    
    
    
    
    
    
    
    
  }


  
  print("all done")
  setwd("/home/master/")
  file.remove(files_to_download2)



}





avg_rasters(flux="NEE_gC_m2", pattern1=seq(1990, 2016, by=1)%>% as.character(), folder="0.5")
avg_rasters(flux="GPP_gC_m2", pattern1=seq(1990, 2016, by=1)%>% as.character(), folder="0.5")
avg_rasters(flux="Reco_gC_m2", pattern1=seq(1990, 2016, by=1)%>% as.character(), folder="0.5")


# years <- seq(1990, 2016, by=1) # TEMPORARY
# 
# for (y in years) {
# 
#   avg_rasters(flux="NEE_gC_m2", pattern1=seq(y, y, by=1)%>% as.character(), folder="0.5")
# 
# }
# 
# for (y in years) {
# 
#   avg_rasters(flux="GPP_gC_m2", pattern1=seq(y, y, by=1)%>% as.character(), folder="0.5")
# 
# }
# 
# for (y in years) {
# 
#   avg_rasters(flux="Reco_gC_m2", pattern1=seq(y, y, by=1)%>% as.character(), folder="0.5")
# 
# }


















### Calcualate annual env summaries


var <- "tmean"
divider <- 1
pattern4 <- seq(1990, 2016) %>% as.character()


avg_maps_env <- function(var, divider, pattern4) {
  
  setwd("/home/master/")
  files_to_download <- grep("*.tif", contents$name, value = TRUE)
  files_to_download2 <- files_to_download[grepl(paste(var), files_to_download)]
  files_to_download2 <- files_to_download2[grepl("predictors_8km", files_to_download2)]
  files_to_download3 <- files_to_download2[!(grepl("predictors_1km", files_to_download2) | grepl("old", files_to_download2))]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2017", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2018", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2019", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2020", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "2021", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "198", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "197", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "196", sep="_"), files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl(paste(var, "195", sep="_"), files_to_download3)) ]
  
  files_to_download3 <- files_to_download3[!(grepl("trend", files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl("_mean_", files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl("_averages_", files_to_download3)) ]
  files_to_download3 <- files_to_download3[!(grepl("/VCF", files_to_download3)) ]
  
  
  map(files_to_download3, function(x) gcs_get_object(x, saveToDisk = x, overwrite = TRUE))
  
  files <- files_to_download3
  
  
  # mean spatial map
  r <- rast(files) # /divider
  rmean <- mean(r, na.rm=TRUE)
  plot(rmean)
  writeRaster(rmean, paste0("/home/master/trends_8km/", var, "_", pattern4[1], "_", pattern4[length(pattern4)],  "_annual_mean.tif"), overwrite=TRUE)
  
  
  # calculate annual means over time
  if (!grepl("Percent", var)) {
    
    r_annuals <- rast()
    ns <- seq(1, 12*length(pattern4), by=12)
    
    for (n in ns) {
      print(n)
      start <- n
      end <- n+11
      
      r_subset <- r[[start:end]]
      print(sources(r_subset))
      
      print(sources(r_subset))
      
      r_subset_mean <- mean(r_subset, na.rm=TRUE)
      
      r_annuals <- c(r_annuals, r_subset_mean)
      
    }
    
    
    
  } else {
    
    r_annuals <- r
    
    
  }
  
  
  
  # loop over years to extract annual means
  rsums <- data.frame()
  for (y in 1: length(pattern4)) {
    
    print(pattern4[y])
    r2 <- r_annuals[[y]]
    rsum_val <- global(r2, "mean", na.rm=TRUE)
    
    bsum_val <- terra::extract(r2, eco, fun="mean", na.rm=TRUE) 
    
    regsum_val <- terra::extract(r2, region, fun="mean", na.rm=TRUE) 
    
    sums_all <- data.frame(rbind(cbind(area="full", mean=unname(rsum_val)), cbind(area=c("Boreal", "Tundra"), mean=bsum_val[,2]), cbind(area=region$region, mean=regsum_val[,2])))
    sums_all$year <- pattern4[y]
    rsums <- rbind(rsums, sums_all)
    
  }
  
  write.csv(rsums, paste0("/home/master/abcfluxv1_modeling/results/", var, "_", pattern4[1], "_", pattern4[length(pattern4)], "_annual_means.csv"))
  
  
  
  
  
  # Time series
  
  rasterdf <-as.data.frame(r_annuals, xy=TRUE)
  trend <- zyp.trend.dataframe(rasterdf, metadata.cols=2, method="yuepilon",
                               conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
  
  trendr <- rast(trend, type="xyz")
  
  writeRaster(trendr, paste0("/home/master/trends_8km/", var, "_", pattern4[1], "_", pattern4[length(pattern4)], "_annual_trend.tif"), overwrite=TRUE)
  
}

# 
# avg_maps_env(var="tmean", divider=1, pattern4=seq(1990, 2016) %>% as.character())
# 
# avg_maps_env(var="ndvi", divider=1, pattern4=seq(1990, 2016) %>% as.character())
# 
# 
# avg_maps_env(var="vpd", divider=1, pattern4=seq(1990, 2016) %>% as.character())
# 
# avg_maps_env(var="snowcover", divider=1, pattern4=seq(1990, 2016) %>% as.character())

# avg_maps_env(var="soilmoist", divider=1, pattern4=seq(1990, 2016) %>% as.character())



#avg_maps_env(var="Percent_NonTree_Vegetation", divider=1, pattern4=seq(1990, 2016) %>% as.character())

avg_maps_env(var="Percent_NonVegetated", divider=1, pattern4=seq(1990, 2016) %>% as.character())


avg_maps_env(var="Percent_TreeCover", divider=1, pattern4=seq(1990, 2016) %>% as.character())

















### calculate budgets and trends from annual files


pattern1 <- "NEE_gC_m2"
pattern2 <- "sum"
pattern3 <- "8km"
pattern4 <- seq(1990, 2016) %>% as.character()
folder <- "0.5"

avg_maps <- function(flux, pattern1, pattern2, pattern3, pattern4, folder) {

  setwd(paste0("/home/master/predictions_8km/raster/", folder))
  files <- list.files(paste0("/home/master/predictions_8km/raster/", folder, "/"), pattern=pattern1)
  files <- files[grepl(pattern2, files)]
  files <- files[grepl(pattern3, files)]
  if(pattern1!="Reco_gC_m2") {
    files <- files[nchar(files)==26]
  } else {
      files <- files[nchar(files)==27]}
  if(pattern1!="Reco_gC_m2") {
    files_years <- substr(files, start=15, stop=18)
  } else {
    files_years <- substr(files, start=16, stop=19)}

  years_select <- files_years  %in% pattern4
  files <- files[years_select]
  print(files)


  # mean spatial map
  r <- rast(files)/100
  rmean <- mean(r, na.rm=TRUE)
  plot(rmean)
  writeRaster(rmean, paste0("/home/master/predictions_8km/raster/", folder, "/", pattern1, "_", pattern3, "_", pattern4[1], "_", pattern4[length(pattern4)], "_", pattern2, ".tif"))

  # mean budget, tundra vs. boreal
  area_raster <- cellSize(rmean, unit="m")
  rsum2 <- rmean*area_raster
  rsum_val <- global(rsum2, "sum", na.rm=TRUE) * 1e-12
  rsum_val


  # loop over years to extract annual budgets
  rsums <- data.frame()
  for (y in 1: length(pattern4)) {

    print(pattern4[y])
    r2 <- r[[y]]
    rsum2 <- r2*area_raster
    rsum_val <- global(rsum2, "sum", na.rm=TRUE) * 1e-12

    bsum_val <- terra::extract(rsum2, eco, fun="sum", na.rm=TRUE) * 1e-12

    regsum_val <- terra::extract(rsum2, region, fun="sum", na.rm=TRUE) * 1e-12

    sums_all <- data.frame(rbind(cbind(area="full", budget=unname(rsum_val)), cbind(area=c("Boreal", "Tundra"), budget=bsum_val[,2]), cbind(area=region$region, budget=regsum_val[,2])))
    sums_all$year <- pattern4[y]
    rsums <- rbind(rsums, sums_all)

  }

  write.csv(rsums, paste0("/home/master/abcfluxv1_modeling/results/", pattern1, "_", pattern3, "_budgets.csv"))





  # Time series

  rasterdf <-as.data.frame(r, xy=TRUE)
  trend <- zyp.trend.dataframe(rasterdf, metadata.cols=2, method="yuepilon",
                               conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)

  trendr <- rast(trend, type="xyz")

  writeRaster(trendr, paste0("/home/master/trends_8km/", folder, "/", pattern1, "_", pattern3, "_", pattern4[1], "_", pattern4[length(pattern4)], "_", pattern2, "_trend.tif"))

}


avg_maps(pattern1="NEE_gC_m2", pattern2="sum", pattern3="8km", pattern4=seq(1990, 2016) %>% as.character(), folder="0.5")


avg_maps(pattern1="GPP_gC_m2", pattern2="sum", pattern3="8km", pattern4=seq(1990, 2016) %>% as.character(), folder="0.5")


avg_maps(pattern1="Reco_gC_m2", pattern2="sum", pattern3="8km", pattern4=seq(1990, 2016) %>% as.character(), folder="0.5")

















# ### calculate pixel-wise environmental trends using zyp package 
# setwd("/home/master/cloud/predictors_8km")
# years <- seq(1990, 2016, 1) %>% as.character()
# vector.is.empty <- function(x) return(length(x) ==0 )
# trends_env <- function(env) {
#   
#   print(env)
#   
#   months <- c("_01", "_02", "_03", "_04", "_05", "_06",
#               "_07", "_08", "_09", "_10", "_11", "_12")
#   
#   files <- list.files("/home/master/cloud/predictors_8km", pattern=env)
#   files <- files[!str_detect(files, "198")] 
#   files <- files[!str_detect(files, "2021")] 
#   files <- files[!str_detect(files, "2022")] 
#   files <- files[!str_detect(files, "trend")] 
#   
#   # # TEMPORARY
#   # files <- files[str_detect(files, "199")] 
#   
#   
#   
#   
#   # annual means
#   
#   
#   
#   # write out table of annual means
#   
#   
#   
#   
#   
#   
#   # trends
#   print(files)
#   
#   
#   for (m in months) {
#     
#     m2 <- sub("0+", "", m)
#     files3 <- files[str_detect(files, m) | str_detect(files, m2)] 
#     
#     r <- rast(files3)
#     df <- as.data.frame(r, xy=TRUE)
#     
#     trend <- zyp.trend.dataframe(df, metadata.cols=2, method="yuepilon",
#                                  conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
#     
#     trendr <- rast(trend, type="xyz")
#     
#     writeRaster(trendr, paste0("/home/master/local_outputs/trends_drivers_8km/", env,  m, "_trend", ".tif"), overwrite=TRUE, datatype="INT4S")
#     
#     
#     print(paste(m, "done"))
#     
#     
#   }
#   
#   
#   # seasons: spring #
#   m <- c("_03_04_05")
#   
#   files3 <- files[str_detect(files, "_03") | str_detect(files, "_04") | str_detect(files, "_05") | str_detect(files, "_3") | str_detect(files, "_4") | str_detect(files, "_5")] 
#   
#   rs <- rast()
#   
#   for (y in years) {
#     
#     print(y)
#     files4 <- files3[str_detect(files3, y) ]
#     r <- rast(files4)
#     r2 <- mean(r)
#     rs <- c(rs, r2)
#     
#   }
#   
#   df <- as.data.frame(rs, xy=TRUE)
#   
#   trend <- zyp.trend.dataframe(df, metadata.cols=2, method="yuepilon",
#                                conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
#   
#   trendr <- rast(trend, type="xyz")
#   
#   writeRaster(trendr, paste0("/home/master/local_outputs/trends_drivers_8km/", env,  m, "_trend", ".tif"), overwrite=TRUE, datatype="INT4S")
#   
#   
#   
#   # seasons: summer #
#   m <- c("_06_07_08")
#   
#   files3 <- files[str_detect(files, "_06") | str_detect(files, "_07") | str_detect(files, "_08") | str_detect(files, "_6") | str_detect(files, "_7") | str_detect(files, "_8")] 
#   
#   rs <- rast()
#   
#   for (y in years) {
#     
#     print(y)
#     files4 <- files3[str_detect(files3, y) ]
#     r <- rast(files4)
#     r2 <- mean(r)
#     rs <- c(rs, r2)
#     
#   }
#   
#   df <- as.data.frame(rs, xy=TRUE)
#   
#   trend <- zyp.trend.dataframe(df, metadata.cols=2, method="yuepilon",
#                                conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
#   
#   trendr <- rast(trend, type="xyz")
#   
#   writeRaster(trendr, paste0("/home/master/local_outputs/trends_drivers_8km/", env,  m, "_trend", ".tif"), overwrite=TRUE, datatype="INT4S")
#   
#   print("spring and summer done")
#   file.remove(list.files("/home/master/temp/", full.names=TRUE))
#   
#   # seasons: autumn #
#   m <- c("_09_10_11")
#   
#   files3 <- files[str_detect(files, "_09") | str_detect(files, "_10") | str_detect(files, "_11") | str_detect(files, "_9") | str_detect(files, "_10") | str_detect(files, "_11")] 
#   
#   rs <- rast()
#   
#   for (y in years) {
#     
#     print(y)
#     files4 <- files3[str_detect(files3, y) ]
#     r <- rast(files4)
#     r2 <- mean(r)
#     rs <- c(rs, r2)
#     
#   }
#   
#   df <- as.data.frame(rs, xy=TRUE)
#   
#   trend <- zyp.trend.dataframe(df, metadata.cols=2, method="yuepilon",
#                                conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
#   
#   trendr <- rast(trend, type="xyz")
#   
#   writeRaster(trendr, paste0("/home/master/local_outputs/trends_drivers_8km/", env,  m, "_trend", ".tif"), overwrite=TRUE, datatype="INT4S")
#   
#   
#   # seasons: winter #
#   m <- c("_01_02_12")
#   
#   files3 <- files[str_detect(files, "_01") | str_detect(files, "_02") | str_detect(files, "_12") | str_detect(files, "_1") | str_detect(files, "_2") | str_detect(files, "_12")] 
#   
#   rs <- rast()
#   
#   years2 <- seq(1989, 2016, 1) %>% as.character()
#   
#   for (y in years2) {
#     
#     print(y)
#     files4 <- files3[str_detect(files3, y) ]
#     files4 <- str_replace_all(files4, paste0(y, "_12"), paste0((as.numeric(y)-1) %>% as.character(), "_12"))
#     r <- rast(files4)
#     r2 <- mean(r)
#     rs <- c(rs, r2)
#     
#   }
#   
#   df <- as.data.frame(rs, xy=TRUE)
#   
#   trend <- zyp.trend.dataframe(df, metadata.cols=2, method="yuepilon",
#                                conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
#   
#   trendr <- rast(trend, type="xyz")
#   
#   writeRaster(trendr, paste0("/home/master/local_outputs/trends_drivers_8km/", env,  m, "_trend", ".tif"), overwrite=TRUE, datatype="INT4S")
#   
#   
#   
#   # seasons: ngs #
#   m <- c("_01_02_03_04_10_11_12")
#   
#   files3 <- files[str_detect(files, "_01") | str_detect(files, "_02") |str_detect(files, "_03") |str_detect(files, "_04") |str_detect(files, "_10") |str_detect(files, "_11") | str_detect(files, "_12") | str_detect(files, "_1") | str_detect(files, "_2") |str_detect(files, "_3") |str_detect(files, "_4")] 
#   
#   rs <- rast()
#   
#   years2 <- seq(1989, 2016, 1) %>% as.character()
#   
#   for (y in years2) {
#     
#     print(y)
#     files4 <- files3[str_detect(files3, y) ]
#     files4 <- str_replace_all(files4, paste0(y, "_12"), paste0((as.numeric(y)-1) %>% as.character(), "_12"))
#     files4 <- str_replace_all(files4, paste0(y, "_11"), paste0((as.numeric(y)-1) %>% as.character(), "_11"))
#     files4 <- str_replace_all(files4, paste0(y, "_10"), paste0((as.numeric(y)-1) %>% as.character(), "_10"))
#     
#     r <- rast(files4)
#     r2 <- mean(r)
#     rs <- c(rs, r2)
#     
#   }
#   
#   df <- as.data.frame(rs, xy=TRUE)
#   
#   trend <- zyp.trend.dataframe(df, metadata.cols=2, method="yuepilon",
#                                conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
#   
#   trendr <- rast(trend, type="xyz")
#   
#   writeRaster(trendr, paste0("/home/master/local_outputs/trends_drivers_8km/", env,  m, "_trend", ".tif"), overwrite=TRUE, datatype="INT4S")
#   
#   
#   # seasons: gs #
#   m <- c("_05_06_07_08_09")
#   
#   
#   files3 <- files[str_detect(files, "_05") | str_detect(files, "_06") | str_detect(files, "_07") | str_detect(files, "_08") | str_detect(files, "_09") | str_detect(files, "_5") | str_detect(files, "_6") | str_detect(files, "_7") | str_detect(files, "_8") | str_detect(files, "_9")] 
#   
#   rs <- rast()
#   
#   for (y in years) {
#     
#     print(y)
#     files4 <- files3[str_detect(files3, y) ]
#     r <- rast(files4)
#     r2 <- mean(r)
#     rs <- c(rs, r2)
#     
#   }
#   
#   df <- as.data.frame(rs, xy=TRUE)
#   
#   trend <- zyp.trend.dataframe(df, metadata.cols=2, method="yuepilon",
#                                conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
#   
#   trendr <- rast(trend, type="xyz")
#   
#   writeRaster(trendr, paste0("/home/master/local_outputs/trends_drivers_8km/", env,  m, "_trend", ".tif"), overwrite=TRUE, datatype="INT4S")
#   
#   
#   # annual #
#   
#   
#   files3 <- files
#   
#   for (y in years) {
#     
#     print(y)
#     files4 <- files3[str_detect(files3, y) ]
#     r <- rast(files4)
#     r2 <- mean(r)
#     rs <- c(rs, r2)
#     
#   }
#   
#   df <- as.data.frame(rs, xy=TRUE)
#   
#   trend <- zyp.trend.dataframe(df, metadata.cols=2, method="yuepilon",
#                                conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
#   
#   trendr <- rast(trend, type="xyz")
#   
#   writeRaster(trendr, paste0("/home/master/local_outputs/trends_drivers_8km/", env,  "_annual","_trend", ".tif"), overwrite=TRUE, datatype="INT4S")
#   
#   
#   rm(trendr);rm(trend);rm(r)
#   gc()
#   file.remove(list.files("/home/master/temp/", full.names=TRUE))
#   
#   
#   
# }
# 
# 
# trends_env("tmean")
# trends_env("ndvi")
# trends_env("srad")
# trends_env("soilmoist")
# trends_env("soiltemp")
# trends_env("vpd")
# 
# 
# 
# 
# 
# # calculate trends of annual conditions
#   
#   
# 
# 
# 
# 
# 
# 
# ### Old not needed because added the folder above
# # 
# # 
# # # 0.025
# # 
# # 
# # 
# # ### 8 km predictions ###
# # 
# # ## In situ data
# # d <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg_allsites.csv")
# # d <- d %>% mutate_if(is.integer, as.numeric)
# # # reverse gpp first, for simplicity
# # d$GPP_gC_m2 <- -d$GPP_gC_m2 
# # 
# # # ### normalize training data (value – min) / (max – min) to be able to scale it back
# # # normalize <- function(x, ...) {
# # #   return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
# # # }
# # # d$NEE_gC_m2 <- normalize(d$NEE_gC_m2, na.rm=TRUE)
# # # d$Reco_gC_m2 <- normalize(d$Reco_gC_m2, na.rm=TRUE)
# # # d$GPP_gC_m2 <- normalize(d$GPP_gC_m2, na.rm=TRUE)
# # 
# # 
# # ### calculate average rasters ###
# # 
# # # remove filename
# # avg_rasters <- function(flux, pattern1) {
# #   
# #   
# #   
# #   setwd("/home/master/")
# #   files_to_download <- grep("*.csv", contents$name, value = TRUE)
# #   files_to_download2 <- files_to_download[grepl(flux, files_to_download)]
# #   files_to_download2 <- files_to_download2[grepl("predictions_8km/csv/0.025/", files_to_download2)]
# #   files_to_download2 <- files_to_download2[!(grepl("train_loocv", files_to_download2)) ]
# #   
# #   
# #   map(files_to_download2, function(x) gcs_get_object(x, saveToDisk = x, overwrite = TRUE))
# #   
# #   
# #   # read to a data frame
# #   setwd("/home/master/predictions_8km/csv/0.025")
# #   files <- list.files(, pattern=flux)
# #   files <- files[!(grepl("train_loocv", files)) ]
# #   df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
# #   
# #   
# #   # rename
# #   columns <- NA
# #   for (f in files) {
# #     column <- c("x", "y", f)
# #     columns <- c(columns, column)
# #     
# #   }
# #   
# #   columns <- columns[2:length(columns)]
# #   names(df) <- columns
# #   
# #   
# #   # normalization values back 
# #   #/ 10000 * (maxval - minval) + minval
# #   # note that the predictions have still been multiplied by 1000
# #   
# #   
# #   minval <- ifelse(flux=="NEE_gC_m2", min(d$NEE_gC_m2, na.rm=TRUE), NA)
# #   minval <- ifelse(flux=="GPP_gC_m2", min(d$GPP_gC_m2, na.rm=TRUE), minval)
# #   minval <- ifelse(flux=="Reco_gC_m2", min(d$Reco_gC_m2, na.rm=TRUE), minval)
# #   
# #   
# #   maxval <- ifelse(flux=="NEE_gC_m2", max(d$NEE_gC_m2, na.rm=TRUE), NA)
# #   maxval <- ifelse(flux=="GPP_gC_m2", max(d$GPP_gC_m2, na.rm=TRUE), maxval)
# #   maxval <- ifelse(flux=="Reco_gC_m2", max(d$Reco_gC_m2, na.rm=TRUE), maxval)
# #   
# #   # normalization values back 
# #   # every third column, do this
# #   cols <- seq(3, length(files)*3, by=3)
# #   df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
# #   
# #   
# #   
# #   
# #   ## temporal splits
# #   
# #   
# #   # months #
# #   months <- c("_01", "_02", "_03", "_04", "_05", "_06",
# #               "_07", "_08", "_09", "_10", "_11", "_12")
# #   
# #   for (m in months) {
# #     # m <- "_01"
# #     selected_columns <- str_detect(columns, m)
# #     
# #     selected_columns[1:2] <- TRUE # add coords
# #     df2 <- df[, selected_columns]
# #     print(colnames(df2)[3:length(colnames(df2))])
# #     
# #     df3 <- df2  # %>% mutate(mean = rowMeans(.[, 3:ncol(df2)])) # this would do annual averages
# #     df4 <- df3[, c(1:2, ncol(df3))]
# #     
# #     df4[, 3] <- unname(df4[, 3])
# #     
# #     rasterdf  <- as.matrix(df4[]) 
# #     r <- rast(rasterdf[], type="xyz")
# #     # plot(r)
# #     crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #     writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", m, "_", pattern1, ".tif"), overwrite=TRUE, datatype="INT4S")
# #     
# #   }
# #   
# #   print("months done")
# #   
# #   # seasons: spring #
# #   months <- c("_03", "_04", "_05")
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S") 
# #   
# #   
# #   
# #   # seasons: summer #
# #   months <- c("_06", "_07", "_08")
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   
# #   
# #   
# #   # seasons: autumn #
# #   months <- c("_09", "_10", "_11")
# #   
# #   
# #   # y <- "2000"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   # y <- "2000"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   
# #   
# #   
# #   print("most of 3-month seasons done")
# #   
# #   
# #   
# #   # Calculate longer seasons (GS)
# #   
# #   # seasons: autumn #
# #   months <- c("_05","_06", "_07", "_08", "_09")
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   print("GS done")
# #   
# #   
# #   
# #   
# #   # annual 
# #   months <- c("_01", "_02", "_03", "_04", "_05", "_06",
# #               "_07", "_08", "_09", "_10", "_11", "_12")
# #   
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   
# #   print("annual done")
# #   
# #   
# #   
# #   # seasons: winter #
# #   # trickier because need to take previous year into account!
# #   months <- c("_01", "_02", "_12")
# #   
# #   
# #   
# #   # list and select files
# #   # files <- list.files("/home/master/predictions_8km/csv/0.025", pattern=paste0(pattern1, collapse="|")) # old way?
# #   files <- list.files("/home/master/predictions_8km/csv/0.025", pattern=pattern1)
# #   files <- files[str_detect(files, flux)]
# #   files <- files[!str_detect(files, "_12")] # remove december from the same year
# #   files2 <- list.files("/home/master/predictions_8km/csv/0.025", pattern=paste0(as.numeric(pattern1)-1, collapse="|"))
# #   files2 <- files2[str_detect(files2, "_12")]
# #   files2 <- files2[str_detect(files2, flux)]
# #   files <- c(files, files2)
# #   files <- files[str_detect(files, flux)]
# #   setwd("/home/master/predictions_8km/csv/0.025")
# #   
# #   # read to a data frame
# #   df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
# #   
# #   
# #   # rename
# #   columns <- NA
# #   for (f in files) {
# #     column <- c("x", "y", f)
# #     columns <- c(columns, column)
# #     
# #   }
# #   
# #   columns <- columns[2:length(columns)]
# #   names(df) <- columns
# #   
# #   
# #   
# #   # normalization values back 
# #   # every third column, do this
# #   cols <- seq(3, length(files)*3, by=3)
# #   df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
# #   
# #   # skip 1982
# #   if (y!="1982") {
# #     
# #     # y <- "1983"
# #     # months 1-2
# #     selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #     
# #     #selected_columns2 <- str_detect(columns, y)
# #     
# #     
# #     #selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE | selected_columnsb==TRUE & selected_columnsb2==TRUE, TRUE, FALSE)
# #     
# #     selected_columns[1:2] <- TRUE # add coords
# #     df2 <- df[, selected_columns]
# #     print(colnames(df2)[3:length(colnames(df2))])
# #     
# #     df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #     
# #     df4 <- df3[, c(1, 2, ncol(df3))]
# #     
# #     rasterdf  <- as.matrix(df4[]) 
# #     r <- rast(rasterdf[], type="xyz")
# #     # plot(r)
# #     crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #     writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #     
# #     
# #     
# #     df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #     
# #     df4 <- df3[, c(1, 2, ncol(df3))]
# #     
# #     rasterdf  <- as.matrix(df4[]) 
# #     r <- rast(rasterdf[], type="xyz")
# #     # plot(r)
# #     crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #     writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #     
# #   }
# #   
# #   
# #   
# #   
# #   # seasons: non-growing season #
# #   # trickier because need to take previous year into account!
# #   months <- c("_01", "_02", "_03", "_04","_10","_11","_12")
# #   
# #   
# #   # list and select files
# #   files <- list.files("/home/master/predictions_8km/csv/0.025", pattern=paste0(pattern1, collapse="|"))
# #   files <- files[!(str_detect(files, "_12") | str_detect(files, "_11") | str_detect(files, "_10"))]# remove december from the same year
# #   files2 <- list.files("/home/master/predictions_8km/csv/0.025", pattern=paste0(as.numeric(pattern1)-1, collapse="|"))
# #   files2 <- files2[str_detect(files2, "_12") | str_detect(files2, "_11") | str_detect(files2, "_10") ]
# #   files <- c(files, files2)
# #   files <- files[str_detect(files, flux)]
# #   setwd("/home/master/predictions_8km/csv/0.025")
# #   
# #   # read to a data frame
# #   df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
# #   
# #   # rename
# #   columns <- NA
# #   for (f in files) {
# #     column <- c("x", "y", f)
# #     columns <- c(columns, column)
# #     
# #   }
# #   
# #   columns <- columns[2:length(columns)]
# #   names(df) <- columns
# #   
# #   
# #   
# #   
# #   # normalization values back 
# #   # every third column, do this
# #   cols <- seq(3, length(files)*3, by=3)
# #   df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
# #   
# #   
# #   for (y in pattern1) {
# #     
# #     # skip 1982
# #     if (y!=1982) {
# #       
# #       # y <- "1983"
# #       # months 1-2
# #       selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #       
# #       #selected_columns2 <- str_detect(columns, y)
# #       
# #       
# #       #selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE | selected_columnsb==TRUE & selected_columnsb2==TRUE, TRUE, FALSE)
# #       
# #       selected_columns[1:2] <- TRUE # add coords
# #       df2 <- df[, selected_columns]
# #       print(colnames(df2)[3:length(colnames(df2))])
# #       
# #       df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #       
# #       df4 <- df3[, c(1, 2, ncol(df3))]
# #       
# #       rasterdf  <- as.matrix(df4[]) 
# #       r <- rast(rasterdf[], type="xyz")
# #       # plot(r)
# #       crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #       writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #       
# #       
# #       df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #       
# #       df4 <- df3[, c(1, 2, ncol(df3))]
# #       
# #       rasterdf  <- as.matrix(df4[]) 
# #       r <- rast(rasterdf[], type="xyz")
# #       # plot(r)
# #       crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #       writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.025/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #       
# #     }
# #     
# #     
# #   }
# #   
# #   
# #   
# #   
# #   print("all done")
# #   setwd("/home/master/")
# #   file.remove(files_to_download2)
# #   
# #   
# #   
# #   
# #   
# # }
# # 
# # 
# # 
# # years <- seq(1990, 2016, by=1) # TEMPORARY
# # 
# # for (y in years) {
# #   
# #   avg_rasters(flux="NEE_gC_m2", pattern1=seq(y, y, by=1)%>% as.character())
# #   
# # }
# # 
# # for (y in years) {
# #   
# #   avg_rasters(flux="GPP_gC_m2", pattern1=seq(y, y, by=1)%>% as.character())
# #   
# # }
# # 
# # for (y in years) {
# #   
# #   avg_rasters(flux="Reco_gC_m2", pattern1=seq(y, y, by=1)%>% as.character())
# #   
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # # 0.975
# # 
# # 
# # 
# # ### 8 km predictions ###
# # 
# # ## In situ data
# # d <- read.csv("/home/master/flux_upscaling_data/results/final/modeldata_avg_allsites.csv")
# # d <- d %>% mutate_if(is.integer, as.numeric)
# # # reverse gpp first, for simplicity
# # d$GPP_gC_m2 <- -d$GPP_gC_m2 
# # 
# # # ### normalize training data (value – min) / (max – min) to be able to scale it back
# # # normalize <- function(x, ...) {
# # #   return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
# # # }
# # # d$NEE_gC_m2 <- normalize(d$NEE_gC_m2, na.rm=TRUE)
# # # d$Reco_gC_m2 <- normalize(d$Reco_gC_m2, na.rm=TRUE)
# # # d$GPP_gC_m2 <- normalize(d$GPP_gC_m2, na.rm=TRUE)
# # 
# # 
# # ### calculate average rasters ###
# # 
# # # remove filename
# # avg_rasters <- function(flux, pattern1) {
# #   
# #   
# #   
# #   setwd("/home/master/")
# #   files_to_download <- grep("*.csv", contents$name, value = TRUE)
# #   files_to_download2 <- files_to_download[grepl(flux, files_to_download)]
# #   files_to_download2 <- files_to_download2[grepl("predictions_8km/csv/0.975/", files_to_download2)]
# #   files_to_download2 <- files_to_download2[!(grepl("train_loocv", files_to_download2)) ]
# #   
# #   
# #   map(files_to_download2, function(x) gcs_get_object(x, saveToDisk = x, overwrite = TRUE))
# #   
# #   
# #   # read to a data frame
# #   setwd("/home/master/predictions_8km/csv/0.975")
# #   files <- list.files(, pattern=flux)
# #   files <- files[!(grepl("train_loocv", files)) ]
# #   df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
# #   
# #   
# #   # rename
# #   columns <- NA
# #   for (f in files) {
# #     column <- c("x", "y", f)
# #     columns <- c(columns, column)
# #     
# #   }
# #   
# #   columns <- columns[2:length(columns)]
# #   names(df) <- columns
# #   
# #   
# #   # normalization values back 
# #   #/ 10000 * (maxval - minval) + minval
# #   # note that the predictions have still been multiplied by 1000
# #   
# #   
# #   minval <- ifelse(flux=="NEE_gC_m2", min(d$NEE_gC_m2, na.rm=TRUE), NA)
# #   minval <- ifelse(flux=="GPP_gC_m2", min(d$GPP_gC_m2, na.rm=TRUE), minval)
# #   minval <- ifelse(flux=="Reco_gC_m2", min(d$Reco_gC_m2, na.rm=TRUE), minval)
# #   
# #   
# #   maxval <- ifelse(flux=="NEE_gC_m2", max(d$NEE_gC_m2, na.rm=TRUE), NA)
# #   maxval <- ifelse(flux=="GPP_gC_m2", max(d$GPP_gC_m2, na.rm=TRUE), maxval)
# #   maxval <- ifelse(flux=="Reco_gC_m2", max(d$Reco_gC_m2, na.rm=TRUE), maxval)
# #   
# #   # normalization values back 
# #   # every third column, do this
# #   cols <- seq(3, length(files)*3, by=3)
# #   df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
# #   
# #   
# #   
# #   
# #   ## temporal splits
# #   
# #   
# #   # months #
# #   months <- c("_01", "_02", "_03", "_04", "_05", "_06",
# #               "_07", "_08", "_09", "_10", "_11", "_12")
# #   
# #   for (m in months) {
# #     # m <- "_01"
# #     selected_columns <- str_detect(columns, m)
# #     
# #     selected_columns[1:2] <- TRUE # add coords
# #     df2 <- df[, selected_columns]
# #     print(colnames(df2)[3:length(colnames(df2))])
# #     
# #     df3 <- df2  # %>% mutate(mean = rowMeans(.[, 3:ncol(df2)])) # this would do annual averages
# #     df4 <- df3[, c(1:2, ncol(df3))]
# #     
# #     df4[, 3] <- unname(df4[, 3])
# #     
# #     rasterdf  <- as.matrix(df4[]) 
# #     r <- rast(rasterdf[], type="xyz")
# #     # plot(r)
# #     crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #     writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", m, "_", pattern1, ".tif"), overwrite=TRUE, datatype="INT4S")
# #     
# #   }
# #   
# #   print("months done")
# #   
# #   # seasons: spring #
# #   months <- c("_03", "_04", "_05")
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S") 
# #   
# #   
# #   
# #   # seasons: summer #
# #   months <- c("_06", "_07", "_08")
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   
# #   
# #   
# #   # seasons: autumn #
# #   months <- c("_09", "_10", "_11")
# #   
# #   
# #   # y <- "2000"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   # y <- "2000"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   
# #   
# #   
# #   print("most of 3-month seasons done")
# #   
# #   
# #   
# #   # Calculate longer seasons (GS)
# #   
# #   # seasons: autumn #
# #   months <- c("_05","_06", "_07", "_08", "_09")
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   print("GS done")
# #   
# #   
# #   
# #   
# #   # annual 
# #   months <- c("_01", "_02", "_03", "_04", "_05", "_06",
# #               "_07", "_08", "_09", "_10", "_11", "_12")
# #   
# #   
# #   
# #   # y <- "1982"
# #   selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #   
# #   selected_columns2 <- str_detect(columns, pattern1)
# #   
# #   selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE, TRUE, FALSE)
# #   
# #   selected_columns[1:2] <- TRUE # add coords
# #   df2 <- df[, selected_columns]
# #   print(colnames(df2)[3:length(colnames(df2))])
# #   
# #   df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #   
# #   df4 <- df3[, c(1, 2, ncol(df3))]
# #   
# #   rasterdf  <- as.matrix(df4[]) 
# #   r <- rast(rasterdf[], type="xyz")
# #   # plot(r)
# #   crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #   writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #   
# #   
# #   
# #   print("annual done")
# #   
# #   
# #   
# #   # seasons: winter #
# #   # trickier because need to take previous year into account!
# #   months <- c("_01", "_02", "_12")
# #   
# #   
# #   
# #   # list and select files
# #   # files <- list.files("/home/master/predictions_8km/csv/0.975", pattern=paste0(pattern1, collapse="|")) # old way?
# #   files <- list.files("/home/master/predictions_8km/csv/0.975", pattern=pattern1)
# #   files <- files[str_detect(files, flux)]
# #   files <- files[!str_detect(files, "_12")] # remove december from the same year
# #   files2 <- list.files("/home/master/predictions_8km/csv/0.975", pattern=paste0(as.numeric(pattern1)-1, collapse="|"))
# #   files2 <- files2[str_detect(files2, "_12")]
# #   files2 <- files2[str_detect(files2, flux)]
# #   files <- c(files, files2)
# #   files <- files[str_detect(files, flux)]
# #   setwd("/home/master/predictions_8km/csv/0.975")
# #   
# #   # read to a data frame
# #   df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
# #   
# #   
# #   # rename
# #   columns <- NA
# #   for (f in files) {
# #     column <- c("x", "y", f)
# #     columns <- c(columns, column)
# #     
# #   }
# #   
# #   columns <- columns[2:length(columns)]
# #   names(df) <- columns
# #   
# #   
# #   
# #   # normalization values back 
# #   # every third column, do this
# #   cols <- seq(3, length(files)*3, by=3)
# #   df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
# #   
# #   # skip 1982
# #   if (y!="1982") {
# #     
# #     # y <- "1983"
# #     # months 1-2
# #     selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #     
# #     #selected_columns2 <- str_detect(columns, y)
# #     
# #     
# #     #selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE | selected_columnsb==TRUE & selected_columnsb2==TRUE, TRUE, FALSE)
# #     
# #     selected_columns[1:2] <- TRUE # add coords
# #     df2 <- df[, selected_columns]
# #     print(colnames(df2)[3:length(colnames(df2))])
# #     
# #     df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #     
# #     df4 <- df3[, c(1, 2, ncol(df3))]
# #     
# #     rasterdf  <- as.matrix(df4[]) 
# #     r <- rast(rasterdf[], type="xyz")
# #     # plot(r)
# #     crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #     writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #     
# #     
# #     
# #     df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #     
# #     df4 <- df3[, c(1, 2, ncol(df3))]
# #     
# #     rasterdf  <- as.matrix(df4[]) 
# #     r <- rast(rasterdf[], type="xyz")
# #     # plot(r)
# #     crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #     writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #     
# #   }
# #   
# #   
# #   
# #   
# #   # seasons: non-growing season #
# #   # trickier because need to take previous year into account!
# #   months <- c("_01", "_02", "_03", "_04","_10","_11","_12")
# #   
# #   
# #   # list and select files
# #   files <- list.files("/home/master/predictions_8km/csv/0.975", pattern=paste0(pattern1, collapse="|"))
# #   files <- files[!(str_detect(files, "_12") | str_detect(files, "_11") | str_detect(files, "_10"))]# remove december from the same year
# #   files2 <- list.files("/home/master/predictions_8km/csv/0.975", pattern=paste0(as.numeric(pattern1)-1, collapse="|"))
# #   files2 <- files2[str_detect(files2, "_12") | str_detect(files2, "_11") | str_detect(files2, "_10") ]
# #   files <- c(files, files2)
# #   files <- files[str_detect(files, flux)]
# #   setwd("/home/master/predictions_8km/csv/0.975")
# #   
# #   # read to a data frame
# #   df <- do.call(cbind,lapply(files,read.csv)) # can do this because the number of rows is same across the predictions!
# #   
# #   # rename
# #   columns <- NA
# #   for (f in files) {
# #     column <- c("x", "y", f)
# #     columns <- c(columns, column)
# #     
# #   }
# #   
# #   columns <- columns[2:length(columns)]
# #   names(df) <- columns
# #   
# #   
# #   
# #   
# #   # normalization values back 
# #   # every third column, do this
# #   cols <- seq(3, length(files)*3, by=3)
# #   df[, cols] <- df[, cols]  / 10000 * (maxval - minval) + minval
# #   
# #   
# #   for (y in pattern1) {
# #     
# #     # skip 1982
# #     if (y!=1982) {
# #       
# #       # y <- "1983"
# #       # months 1-2
# #       selected_columns <- str_detect(columns, paste(months, collapse = "|")) 
# #       
# #       #selected_columns2 <- str_detect(columns, y)
# #       
# #       
# #       #selected_columns <- ifelse(selected_columns==TRUE & selected_columns2==TRUE | selected_columnsb==TRUE & selected_columnsb2==TRUE, TRUE, FALSE)
# #       
# #       selected_columns[1:2] <- TRUE # add coords
# #       df2 <- df[, selected_columns]
# #       print(colnames(df2)[3:length(colnames(df2))])
# #       
# #       df3 <- df2  %>% mutate(sum = rowSums(.[, 3:ncol(df2)]))
# #       
# #       df4 <- df3[, c(1, 2, ncol(df3))]
# #       
# #       rasterdf  <- as.matrix(df4[]) 
# #       r <- rast(rasterdf[], type="xyz")
# #       # plot(r)
# #       crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #       writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_sum.tif"), overwrite=TRUE, datatype="INT4S")
# #       
# #       
# #       df3 <- df2  %>% mutate(sum = rowMeans(.[, 3:ncol(df2)]))
# #       
# #       df4 <- df3[, c(1, 2, ncol(df3))]
# #       
# #       rasterdf  <- as.matrix(df4[]) 
# #       r <- rast(rasterdf[], type="xyz")
# #       # plot(r)
# #       crs(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# #       writeRaster(r*100, paste0("/home/master/predictions_8km/raster/0.975/", flux, "_8km", stri_paste(months, collapse=''), "_", pattern1, "_mean.tif"), overwrite=TRUE, datatype="INT4S")
# #       
# #     }
# #     
# #     
# #   }
# #   
# #   
# #   
# #   
# #   print("all done")
# #   setwd("/home/master/")
# #   file.remove(files_to_download2)
# #   
# #   
# #   
# #   
# #   
# # }
# # 
# # 
# # 
# # years <- seq(1990, 2016, by=1) # TEMPORARY
# # 
# # for (y in years) {
# #   
# #   avg_rasters(flux="NEE_gC_m2", pattern1=seq(y, y, by=1)%>% as.character())
# #   
# # }
# # 
# # for (y in years) {
# #   
# #   avg_rasters(flux="GPP_gC_m2", pattern1=seq(y, y, by=1)%>% as.character())
# #   
# # }
# # 
# # for (y in years) {
# #   
# #   avg_rasters(flux="Reco_gC_m2", pattern1=seq(y, y, by=1)%>% as.character())
# #   
# # }
# # 
# # 
# # 
