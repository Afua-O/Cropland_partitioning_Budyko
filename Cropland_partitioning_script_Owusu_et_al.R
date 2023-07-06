##R script for Owusu et al
#A Framework for Disaggregating Remote-Sensing Cropland into Rainfed and 
#Irrigated Classes at Continental Scale 

library(dplyr)
library(terra)
library(rasterVis)
library(rts)
library(lubridate)
library(sf)
library(gtools)
library(Rcpp)
library(exactextractr)

#_______________________________________________________________________________
#Set working directory
rm(list=ls()) #remove existing files
setwd("/home/tmp_disk/")
print(getwd())

#For creating annual data stacks
date_month<- seq(as.Date("2019/01/01"), as.Date("2019/12/31"), by="month")
date_year<- seq(as.Date("2019/01/01"), as.Date("2019/12/31"), by="year")

#_______________________________________________________________________________
#Import climate variables:  P, ET, ETo

#P
path = "/home/tmp_disk/P_Africa/Mswep"
P_list <- list.files(path=path, pattern='tif$', full.names=TRUE)
P_stack <- rast(P_list)

#ET 
path = "/home/tmp_disk/ET_Africa/Ssebop/"
ET_list <- list.files(path=path, pattern='tif$', full.names=TRUE)
ET_stack <- rast(ET_list)

#ETo
path = "/home/tmp_disk/ETo_Africa/WaporEto"
ETo_list <- list.files(path=path, pattern='tif$', full.names=TRUE)
ETo_stack <- rast(ETo_list)

#_______________________________________________________________________________
#Import probability cropland map of highest performing cropland masks and derive HCCM

path = "/home/tmp_disk/LU_Africa/Crop/"
LU_list <- list.files(path=path, pattern='tif$', full.names=TRUE)
LU_stack <- rast(LU_list)
LU_stack[LU_stack>=60]<- 1
LU_stack[LU_stack!=1]<- 0
names(LU_stack) <- date_year

#_______________________________________________________________________________
#Aggregate climate inputs to annual timestep

#Precipitation 
P_ts <- rts(P_stack, date_month)
P_ts <- apply.yearly(P_ts, sum) #raster timeseries
P_year1 <-  P_ts@raster #raster stack
names(P_year1) <- date_year

#Actual ET (Using Ssebop Ssebop annual data for Africa)
ET_ts <- rts(ET_stack, date_month)
ET_ts <- apply.yearly(ET_ts, sum) #raster timeseries
ET_year <-  ET_ts@raster #raster stack
ET_year <- ET_stack
names(ET_year) <- date_year

#Potential ET
ETo_ts <- rts(ETo_stack, date_month)
ETo_ts <- apply.yearly(ETo_ts, sum) #raster timeseries
ETo_year <-  ETo_ts@raster #raster stack
names(ETo_year) <- date_year

#_______________________________________________________________________________
#Bias correction of P 

pCF = 1/0.95
P_year <-  P_year1 * pCF 

#_______________________________________________________________________________
#Resampling to 300m  using landcover product:
P_year<- resample(P_year, LU_stack)
P_year

ET_year <- resample(ET_year, LU_stack)
ET_year

ETo_year <- resample(ETo_year, LU_stack)
ETo_year

#_______________________________________________________________________________
# #Import LULC with WA classes and extract natural grasslands

path = "/home/tmp_disk/LU_Africa/Grass/"
Grss_list <- list.files(path=path, pattern='tif$', full.names=TRUE)
Grss <- rast(Grss_list)
Grss[Grss>=80]<- 1
Grss[Grss!=1]<- 0
names(Grss) <- date_year

#_______________________________________________________________________________
#Error correction of ET

P_list <- as.list(P_year)
ET_list <- as.list(ET_year)
Grss_list <- as.list(Grss)

#-------------------------------------------
#Extract ET for grassland for each month

for (i in 1:length(ET_list)){
  
  x <- rast(ET_list[i])
  #print(x)
  y <- rast(Grss_list[i])
  #print(y)
  ET_Grss <- selectRange(x,y)
  writeRaster(ET_Grss, paste("/home/tmp_disk/Output/output_ET_Grss/", "year_", i, ".tif", sep = ""), overwrite=TRUE)
} 
ETGrss_list <- list.files("/home/tmp_disk/Output/output_ET_Grss/", pattern =  ".tif", full.names = T)
ETGrss_stack <- rast(mixedsort(ETGrss_list))

#-------------------------------------------
#Extract P for grassland for each month

for (i in 1:length(P_list)){
  
  x <- rast(P_list[i])
  #print(x)
  y <- rast(Grss_list[i])
  #print(y)
  P_Grss <- selectRange(x,y)
  writeRaster(P_Grss, paste("/home/tmp_disk/Output/output_P_Grss/", "year_", i, ".tif", sep = ""), overwrite=TRUE)
} 
PGrss_list <- list.files("/home/tmp_disk/Output/output_P_Grss/", pattern =  ".tif", full.names = T)
PGrss_stack <- rast(mixedsort(PGrss_list))

#-------------------------------------------
Diff <-   ETGrss_stack - PGrss_stack
plot(Diff, range = c(-1000, 1000), col = rev(hcl.colors(8, "Spectral", rev=FALSE)))
Diff[Diff<=0]= NA 

#Extract mean Diff per basin
Africa <- st_read("/home/tmp_disk/shapefile/Africa_subbasins.shp")
Africa <- Africa[, -(2:15)] 
Err_mean <- exact_extract(Diff, Africa, "mean")

#Attach the mean error value to the attributes of the shapefile
Africa1 <-Africa%>%mutate(Err_mean)

#Rasterize polygons
x2019 <- rasterize(Africa1, ETo_year, field="Err_mean", fun = max)

CF <- c(x2019)
CF[is.nan(CF)] <-0
plot(is.nan(CF)) 
#_______________________________________________________________________________
#apply error correction to ET_stack

zero <- ET_year - ET_year 
ET_yearCF <- max(zero, ET_year - CF) 

#_______________________________________________________________________________
#Calculate green and blue ET using Budyko model

aridi <- ETo_year/P_year
EPratio <- sqrt(aridi*(tanh(1/aridi)* (1-exp(-aridi))))
ETg <- min(1 * EPratio * P_year, ET_yearCF) #no factor applied (factor = 1)
ETb <- ET_yearCF - ETg

#_______________________________________________________________________________
# Create a list of annual ETb and crop LU rasters

ETb_list <- as.list(ETb)
names(LU_stack) <- date_year
plot(LU_stack, col = hcl.colors(14,"Earth", rev = FALSE))
LU_list <- as.list(LU_stack)

#_______________________________________________________________________________
#Calculate subbasin threshold of ETb for irrigated croplands

for (i in 1:length(ETb_list)){
  
  x <- rast(ETb_list[i])
  #print(x)
  y <- rast(Grss_list[i])
  #print(y)
  ETb_Grss <- selectRange(x,y)
  writeRaster(ETb_Grss, paste("/home/tmp_disk/Output/output_ETb_Grss/", "year_", i, ".tif", sep = ""), overwrite=TRUE)
} 
ETbGrss_list <- list.files("/home/tmp_disk/Output/output_ETb_Grss/", pattern =  ".tif", full.names = T)
ETbGrss_stack <- rast(mixedsort(ETbGrss_list))
names(ETbGrss_stack)<-date_year

#-------------------------------------------------------------------------------
#Find average threshold value
Thrshld_mean<- exact_extract(ETbGrss_stack, Africa, "mean")

#Attach the mean error value to the attributes of the shapefile
Africa
Africa2 <-Africa%>%mutate(Thrshld_mean)

#Rasterize polygons
y2019 <- rasterize(Africa2, ETo_year, field="Thrshld_mean", fun = max)
Thrshld <- c(y2019)
Thrshld[is.nan(Thrshld)] <-0
plot(is.nan(Thrshld))

#_______________________________________________________________________________
#Extract ETb for cropland areas in LU maps (ie: where LU ==1)

for (i in 1:length(ETb_list)){
  
  x <- rast(ETb_list[i])
  y <- rast(LU_list[i])
  ETb_LU <- selectRange(x,y)
  writeRaster(ETb_LU, paste("/home/tmp_disk/Output/output_ETb_cropLU/", "year_", i, ".tif", sep = ""), overwrite=TRUE)
}
ETbLU_list <- list.files("/home/tmp_disk/Output/output_ETb_cropLU/", pattern =  ".tif", full.names = T)
ETbLU_stack <- rast(mixedsort(ETbLU_list))
plot(ETbLU_stack,range = c(0, 2000), col = hcl.colors(8, "Oslo", rev = FALSE))

#_______________________________________________________________________________
#Identification of croplands with ETb > threshold
ETb_irr <- ETbLU_stack >= (max(Thrshld,0.01))
plot(ETb_irr,range = c(0, 2000), col = hcl.colors(100, "BluYl", rev = TRUE))
freq(ETb_irr, value=NULL, bylayer=TRUE, usenames=TRUE)

setwd("/home/tmp_disk/Output/output_irri_Budyko_only/")

#create list of ETb_irr
path = getwd()
ETbIrr_list <- list.files(path = path, pattern='tif$', full.names=TRUE)
ETbIrr<- rast(ETbIrr_list)

#_______________________________________________________________________________
#Identify cropland areas of high ETb (formal irrigation) vs low ETB (supplementary irrigation)
formal_irri<-ETbLU_stack
formal_irri[formal_irri<100] <-0
formal_irri[formal_irri>=100] <-1

plot(formal_irri)
setwd("/home/tmp_disk/Output/output_formal_irrigation/")
writeRaster(formal_irri,paste0(names(formal_irri),".tif"), overwrite=TRUE)

##_______________________________________________________________________________
# Extract stats

#Import Africa shapefiles (Basin, wholw continent, Countries)
AfricaB <- st_read("/home/tmp_disk/shapefile/Africa_basins.shp")
Africa_whole <- st_read("/home/tmp_disk/shapefile/Africa_Madagascar.shp") # Continental Africa
Africa_C <- st_read("/home/tmp_disk/shapefile/Africa_countries.shp") # Africa countries

#Stats for cropland vs non-cropland
#number of cropland pixels- continental
exact_extract(LU_stack, st_as_sf(Africa_whole), c('count', 'frac'), append_cols = 'Id') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n)) 

#number of cropland pixels- country
exact_extract(LU_stack, st_as_sf(Africa_C), c('count', 'frac'), append_cols = 'ADM0_NAME') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n))

#number of cropland pixels- basin
exact_extract(LU_stack, st_as_sf(AfricaB), c('count', 'frac'), append_cols = 'HYBAS_ID') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n))

#-----
#Stats for Budyko irrigated vs non-irrigated
#number of pixels- continental
exact_extract(ETbIrr, st_as_sf(Africa_whole), c('count', 'frac'), append_cols = 'Id') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n)) 

#number of pixels- country
exact_extract(ETbIrr, st_as_sf(Africa_C), c('count', 'frac'), append_cols = 'ADM0_NAME') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n))

#number of pixels- basin
exact_extract(ETbIrr, st_as_sf(AfricaB), c('count', 'frac'), append_cols = 'HYBAS_ID') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n))

#---------------
#Stats for formal irrigation irrigated vs supplementary/informal irrigation
#number of pixels- continental
exact_extract(formal_irri, st_as_sf(Africa_whole), c('count', 'frac'), append_cols = 'Id') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n)) 

#number of pixels- country
exact_extract(formal_irri, st_as_sf(Africa_C), c('count', 'frac'), append_cols = 'ADM0_NAME') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n))

#number of pixels- basin
exact_extract(formal_irri, st_as_sf(AfricaB), c('count', 'frac'), append_cols = 'HYBAS_ID') %>%
  mutate(across(starts_with('frac'), function(x) x * count))  %>%
  rename_with(function(n) sub('frac', 'freq', n))
