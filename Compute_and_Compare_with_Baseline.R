#install.packages(c("raster", "rgdal", "gdalUtils", "maptools", "sp", "mapview"))

library(raster)
library(rgdal)
library(gdalUtils)
library(maptools)
library(sp)
library(mapview)
#library(parallel)
#library(doParallel)


#### Workingdirectory ####
#### *** Set your working dir ***
setwd("")


#### Inputdata ####

## loads the point shape for extracting the rastervalues 
shp<-readOGR(paste0(getwd(),"/Points_khorezm.shp")) ##point shapefile

#NDVI Data
#Data for Baseline
stack_ndvi <- stack(paste0(getwd(),"/NDVI_Khorezm.grd"))

#Data with current Images (Images to compare)
current_images_NDVI <- brick(paste0(getwd(),"/current_NDVI_Khorezm.grd"))

#ET Data
#Data for Baseline
stack_ET <- stack(paste0(getwd(),"/ET_Khorezm.grd"))
current_images_ET <- brick(paste0(getwd(),"/current_ET_Khorezm.grd"))




### Preprocess Data ####

Stack_Resample <- function(images, res = 5000){
  
  #images: list of Images to be stacked
  #res   : Resolution for resampling raster cells; Set to 5 km as required 
  #crs   : Crs-string of target crs (Sinusoidal- Proj: "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")
  
  
  
  ## stacking all layers
  ndvi_brick<- brick(images) 
  
  #Extract CRS
  crs1 <- ndvi_stack@crs
  
  ## creating a grid as Mask for resampling -> resampled imageresolution is supposed to be 5x5 km 
  grid_Raster <- raster(ext=extent(ndvi_brick), resolution=res, crs= crs1) # the crs can be extracted from the stack or one single image
  writeRaster(grid_Raster, filename = paste0(getwd(),"/grid_Raster.tif"), format="GTiff")
  
  ## resampling and saving the stack 5x5 km 
  ndvi_brick <- resample(x=ndvi_stack, y=grid_Raster)
  #saving brick
  writeRaster(ndvi_brick, filename='D:UNI/Erstellte_Daten/ndvi_brick_khorezm.tif', format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW")) #saving brick on local drive
  
  return(ndvi_brick)
}



#### Defining the Functions to process data Process and finally combining them into one Function ####

### Convert from Raster to point df ###

Extract_to_points <- function(brick, shp = shp){
  
  #brick : Brick for baseline extraction
  #shape : Spatialpoints df of the study area; Default: shp = Khorezm region
  
  
  #Extracts NDVI-Values at every point location from the providet shp for each layer in the stack 
  point_value<- extract(brick, shp, fun=mean, df=T, na.rm=T)
  
  #saving dataframe
  saveRDS(point_value, file="point_value.Rda")  
  
  return(point_value)
}


### Calc Baseline from reference period ###

Baseline <- function(point_value, n, k, var){
  
  #point_value              : Data.frame containing Spatial points with values corrosponding to Rastercellvalues and Layers 
  #n                        : number of columns in dataframe having images information plus 1 for each year of baseperiod (2001-2018)  
  #k                        : first day of data in each year
  #var                      : Defines the variable witch is calculated mean, sd or quantile 
  #i & j                    : Column number of start of data in dataframe 
  
  
  
  #Variabels for-loop
  i = 2
  j = 2
  
  
  #create empty Data.frame to store results 
  rows <- nrow(point_value)
  newdataframe<- data.frame(c(1:rows))
  
  # change calculation (apply fun) & collumn name if var == sd
  msd <- "median"
  x_string <- "Median_Day_"
  
  if(var == "m"){
    msd <- "median"
    x_string <- "Median_Day_"
  }else{
    msd <- "sd"
    x_string <- "SD_Day_"
  }
  
  for(i in 2:n) {
    #Define Index for Median Calc to extract values from point df.Parameter def.: i ==  829 == all Images (Baseline) by=46 == 1 Year   
    x <- seq(i, 460, by=46)   #there are 829 images in total from 2001-2018 baseline, if need to chnage baseline then number can be changed
    #Extraxt Median values based on x (Index Value )
    dat<- point_value[x]
    #row.names(dat)
    #colnames(dat)
    d<- apply(dat, 1, msd)    #  Use 1 to perform arithmatic on dataframe rows and 2 for columns
    d <- as.data.frame(d)
    
    index <- j      # identifying column where to store the results
    day <- k
    final <- cbind(newdataframe, d)
    colnames(final)[index] <- paste(x_string, " ", day)
    newdataframe <- final
    #View(final)
    #View(dat)
    j <- j+1
    k <- 8+k
  }
  
  return(final)
}


### back to rasters ###

ResultRaster <- function(shp = shp, baseline_median, baseline_sd, stack){
  
  #converts the spatial data.frame back to a Raster Dataset
  
  #shp                      : Studyarea shp with xy-coords
  #baseline_median          : Baseline median from function Baseline
  #baseline_sd              : Baseline sd from function Baseline 
  #stack                    : NDVI_Stack to get original crs 
  
  
  #Bind the median and sd data.frame together with coordinates from the .shp 
  result_table <- cbind(shp$xcoord, shp$ycoord, baseline_median[,2:47], baseline_sd[,2:47])
  
  #Create a Raster with median Values
  dfr1 <- rasterFromXYZ(result_table[,1:48])   #saves all data columns into separate layers as a brick
  
  #Create a Raster with sd Values
  dfr2 <- rasterFromXYZ(result_table[,c(1:2,49:94)])
  
  #dfr <- rasterFromXYZ(make_shape@data[,1:48])   #saves all data columns into separate layers as a brick
  crs(dfr1) <- crs(stack) #defining projetion 
  crs(dfr2) <- crs(stack) #defining projetion 
  
  #calc upper Limit (median + sd)
  sd_upper_limb <- (dfr1)+(dfr2)
  
  #calc lower Limit  (median - sd)
  sd_lower_limb <- (dfr1)-(dfr2)
  
  
  
  #Store results in a list object and return the list (multiple retuns not allowed)
  result_list <- list(Median_df = dfr1, Sd_df = dfr2, sd_upper_limb = sd_upper_limb, sd_lower_limb = sd_lower_limb)
  
  return(result_list)
}



### Compare Baseline with current Image ###
#set desired image of year
#Timestep considered for drought detection (one out of 46 images in one year)
Comp_with_Baseline  <- function(imageofyear, current_images, sd_upper, sd_lower, fun){
  
  #image of year          : Image of the year out of total 46 whose drought map is to be produced
  #current_images         : list of current images 
  #sd_upper               : Raster containing the upper sd values from Baseline function
  #sd_lower               : Raster containing the lower sd values from Baseline function
  
  
  #current_ndvi_brick<- current_images[1:imageofyear]
  #ndvi_current     <- ndvi_brick[[829:874]]   #images of last year whose drought to be predicted i.e. 2019
  
  current_ndvi <- current_images[[imageofyear]]
  sd_upper_baseline  <- sd_upper[[imageofyear]]
  sd_lower_baseline  <- sd_lower[[imageofyear]]
  
  # print(paste0(current_ndvi,"  current image"))
  # print(paste0(sd_upper_baseline, " sd_upper"))
  # print(paste0(sd_lower_baseline, "sd_lower"))
  
  #Stacking all layers included in comparison
  comp_stack  <- stack(current_ndvi, sd_upper_baseline, sd_lower_baseline)
  
  #print(paste0(comp_stack, "Final stack"))
  
  #compare current Timestep with Upper and Lower limits of normal conditions
  
  conditional_statement_drought <- function(current_NDVI, sd_upper_baseline, sd_lower_baseline){
    
    ifelse(current_NDVI > sd_upper_baseline, 1, ifelse( current_NDVI > sd_lower_baseline, 0, -1))
  }
  
  #Variable name needs to be changed 
  Drought_Index_Map <- overlay(comp_stack, fun = conditional_statement_drought)
  
  return(Drought_Index_Map)
}




## outsourced calculations from "WL_NDVI"-Function, because of runtime issues (took to long if calculated in "WL_NDVI"-Function)
# calc Baseline for median and sd
ndvi_points <- Extract_to_points(stack_ndvi, shp)
base_meadian <- Baseline(point_value = ndvi_points, n = 47, k = 1, var = "m")
base_sd      <- Baseline(point_value = ndvi_points, 47, k = 1, var = "sd")


## Final WL Function (combining the defined functions from above)
WL_NDVI <- function(shp, ndvi_stack,base_meadian, base_sd, current_images, current_image_year){
  
  #shp                 : Shp of the desiered study area 
  #ndvi_stack          : Ndvi_stak for Baseline calculation
  #current_images      : Stack with available current images e.g all images of 2019 (CompwithBaseline function)
  #current_image_year  : Desired image out of the current year stack witch is compared against Baseline  
  
  if(current_image_year < 9 || current_image_year >= 39)
  {
    stop("The selcted image is out of bounds")
  }else
  {
    next
  }
  
  
  ## we are working with a subseted test Dataset, thats why we dont need to run the next couple of lines
  #test Data is already resampled
  #resampled_NDVI <- Stack_Resample(stack_ndvi, res = 5000)
  #ndvi_points <- Extract_to_points(stack_ndvi, shp)
  
  ## to speed up the calculation time on our netbooks, we run the next lines apriori 
  # #calc Baseline for median and sd
  # base_meadian <- Baseline(point_value = ndvi_points, n = 47, k = 1, var = "m")
  # base_sd      <- Baseline(point_value = ndvi_points, 47, k = 1, var = "sd")
  
  ## in case you wish to check the standard deviation
  #print(base_sd)
  
  
  #create Rasterdata from pointdata 
  ResultRasterTest <- ResultRaster(shp,
                                   baseline_median = base_meadian,
                                   baseline_sd = base_sd,
                                   stack =  stack_ndvi)
  
  #Creating variabels from Result list
  sd_upper <- ResultRasterTest$sd_upper_limb
  sd_lower <- ResultRasterTest$sd_lower_limb 
  
  #print(sd_lower)
  
  #Comparisn between Baseline and current Timestep
  TestCompwithBase <- Comp_with_Baseline(imageofyear = current_image_year ,current_images = current_images,
                                         sd_upper = sd_upper,
                                         sd_lower = sd_lower)
  
  ## in case you want the check the comarision between the Base line and the thresholds 
  #plot(TestCompwithBase)
  #get_pixel_value <- click(TestCompwithBase, cell=TRUE)
  
  return(TestCompwithBase)
  
}



## End of defining the Functions 
## Feed the Function ####


Wl_stack <- brick()
x <- 39
i <- 0

# 
for (i in 9:x) {
  
  # x       : Betrachtete Zeitschritt
  # x < 9   : Abfrage ob image of year < 06.03 <- Beginn der Beobachtungsperiode
  # x >= 39 : Abfrage ob image of year >= 01.11 <- Ende der Beobachtungsperiode
  
  if(x < 9 || x > 39){
    
    break 
    
  }
  else
  {
    
    store <- WL_NDVI(shp = shp,
                     ndvi_stack = stack_ndvi,
                     base_meadian = base_meadian,
                     base_sd =base_sd, 
                     current_images = current_images_NDVI, 
                     current_image_year = i)
    
    Wl_stack_1 <- addLayer(Wl_stack_1, store)
  }
  
}

for (i in 9:x) {
  
  # x       : Betrachtete Zeitschritt
  # x < 9   : Abfrage ob image of year < 06.03 <- Beginn der Beobachtungsperiode
  # x >= 39 : Abfrage ob image of year >= 01.11 <- Ende der Beobachtungsperiode
  
  if(x < 9 || x > 39){
    
    break 
    
  }
  else
  {
    
    store <- WL_NDVI(shp = shp,
                     ndvi_stack = stack_ET,
                     base_meadian = base_meadian,
                     base_sd =base_sd, 
                     current_images = current_images_ET, 
                     current_image_year = i)
    
    Wl_stack_2 <- addLayer(Wl_stack_2, store)
  }
  
}


Wl_stack_NDVi <- Wl_stack_1
Wl_stack_ET <- Wl_stack_2

writeRaster(Wl_stack_NDVi, "Wl_stack_NDVi.grd", format = "raster", overwrite = TRUE)
writeRaster(Wl_stack_ET, "Wl_stack_NDVi.grd", format = "raster", overwrite = TRUE)


