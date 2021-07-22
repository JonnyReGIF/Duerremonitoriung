
########################################
#### R U L S E T   F U N C T I O N S ####
########################################

## Set working dir
## *** set your working dir ***  
setwd()

## Inputdata 
Wl_stack_NDVi <- stack(paste0(getwd(), "/Wl_stack_NDVi.grd"))
Wl_stack_ET   <- stack(paste0(getwd(), "/Wl_stack_ET.grd"))



## analyses and classification of the current drought situation 
## Function to rename the layer names into their specific dates
dates <- function(Wl_stack_ET) {
  a <- seq(as.Date("06-03", format="%d-%m"), as.Date("01-11", format="%d-%m"), by="days")
  b <- a[seq(1, length(a), 8)]
  c <- b[1:nlayers(Wl_stack_ET)]
  c <- paste("Datum ", format(c, format="%d-%B"))
  return(c)}

##Functions to analyse Drougthcondition



## Function to calculate the DSI

Fun_DSI <- function(Wl_stack_NDVi, Wl_stack_ET, imageofyear){
  
  g <- brick()
  
  for ( i in 1:nlayers(Wl_stack_NDVi)) {
    h <- Wl_stack_NDVi[[i]] + Wl_stack_ET[[i]] 
    
    g <- addLayer(g, h)
    
  }
  #names(g) <- dates(Wl_stack_NDVi)
  names(g) <- dates(Wl_stack_ET)
  return(g) 
}



DSI <- Fun_DSI(Wl_stack_NDVi, Wl_stack_ET, 20)

plot(DSI)


## Function to calculate the drought situation for each timestep 
# for each pixel in each timestep the drought situation is assigned with "0" = no drought OR "1" = drought
# therefore the DSI is used 

Fun_Drought_Situation <- function(DSI){
  g <- brick()
  
  for (i in 1:nlayers(DSI)) {
    
    h <- overlay(DSI[[i]], fun =
                   function(DSI){
                     ifelse(DSI >= 0, 0, 1)
                   })
    
    g <- addLayer(g, h)
  }
  names(g) <- dates(DSI)
  return(g)
}



Drought_Situation <- Fun_Drought_Situation(DSI)

plot(Drought_Situation)




## Function to sum up the Duration of the drought 
# To restart after the end of a drought we mask the respective layer. Due to that the pixel of the summed up Layer will become Zero if the the current pixel is Zero (e.g. no Drought)  

Fun_Duration_Sum <- function(Drought_Situation) {
  g <- brick() ## creats an emty brick to save the summed up data in 
  g <- addLayer(g, Drought_Situation[[1]]) ## adds the first layer of the Current_Drought_Situation-Layer to the new brick (because there is no data from timesteps before there is no need to sum up this information)
  
  for (i in 2:nlayers(Drought_Situation)) {   ## the for loop fills the new brick with layers containing the summed up data for each time step 
    m <- mask(g[[i-1]], Drought_Situation[[i]], maskvalue = 0, updatevalue = 0) ## the mask command resets the those Pixels who turn into zero. If there is no drought detected in a pixel of the 4th timestep we use that to reset the same pixel of the summed-up-Layer from timestep 3 to Zero. All other pixels keep there summed up vlaues.
    j <- Drought_Situation[[i]]+m ## now we can add the manipulated masked layer to the the current drought situation layer, getting a new layer where every pixel where currently is no Drought contains a Zero. all others are summed up
    g <- addLayer(g, j) ## adding the new layer to the brick
  }
  names(g) <- dates(Drought_Situation)
  return(g) ## returns the filled brick
}


Duration_Sum <- Fun_Duration_Sum(Drought_Situation)

plot(Duration_Sum)




### DSI-Sum Function
# 


Fun_DSI_Sum <- function(DSI){
  values(DSI)[values(DSI) >= 0] = 0 # all values that are greater or equal than 0 will be replace with Zero
  g <- brick()
  g <- addLayer(g, DSI[[1]])
  
  
  for (i in 2:nlayers(DSI)) {
    m <- mask(g[[i-1]], DSI[[i]], maskvalue = 0, updatevalue = 0)  # to stop calculation when a pixel will turn into Zero, we manipulate layer from last timestep and imprint a the zeros from current timestep
    j <- DSI[[i]]+m
    g <- addLayer(g, j)
  }
  names(g) <- dates(DSI)
  return(g)
}


DSI_Sum <- Fun_DSI_Sum(DSI)
plot(DSI_Sum)


## Function to calculate Severity
# therefor the Sum of DSI and Duration is needed
# formular was provided  

Fun_Severity <- function(DSI_Sum, Duration_Sum) {
  g <- brick()
  
  for (i in 1:nlayers(DSI_Sum)) {
    h <- DSI_Sum[[i]]
    d <- Duration_Sum[[i]]
    
    j <- overlay(h, d, fun = function(h, d){
      ifelse(d <= 0, 0, abs(h)/d*1.5)})
    g <- addLayer(g, j)
    
  }
  g[is.nan(g)] <- 0 
  names(g) <- dates(DSI_Sum)
  return(g)
}



Severity <- Fun_Severity(DSI_Sum, Duration_Sum)
plot(Severity[[2]])



## Function for the classification of drought

## classification and  encoding of the different drought situation classes
# Water surplus = 1
# No drought = 2
# initial mild drought = 3
# initial severe drought = 4
# mid-term mild drought = 5
# mid-term severe drought = 6
# long term mild drought = 7
# long severe drought = 8


## 
Fun_Class_Drought <- function(Wl_stack_ET, Duration_Sum, Severity){
  g <- brick()
  
  for (i in 1:nlayers(Wl_stack_ET)) {
    wl <- Wl_stack_ET[[i]]
    d <- Duration_Sum[[i]]
    s <- Severity[[i]]
    
    j <- overlay(wl, d, s, fun = 
                   function(wl, d, s){
                     ifelse(wl>=1, 1, 
                            ifelse(d == 0, 2, 
                                   ifelse(d >= 1 & d <= 2 & s < 1, 3,
                                          ifelse(d >= 1 & d <= 2 & s >= 1, 4,
                                                 ifelse(d >= 2 & d <= 3 & s < 1, 5, 
                                                        ifelse(d >= 2 & d <= 3 & s >= 1, 6,
                                                               ifelse(d > 3 & s < 1, 7, 
                                                                      ifelse(d > 3 & s >= 1, 8, 0))))))))}) 
    g <- addLayer(g, j)
  }
  names(g) <- dates(Wl_stack_ET)
  return(g)
  
}


Classification <- Fun_Class_Drought(Wl_stack_ET, Duration_Sum, Severity) 
#plot(Classification[[12]], main = "Classification of the drought situation") 
writeRaster(Classification, "Classes.grd", format = "raster", overwrite = TRUE)

## reproject raster to utm zone 40 
Classification_UTM <- projectRaster(Classification, crs = "EPSG:32640", method = "ngb")
writeRaster(Classification_UTM, "Classes_UTM.grd", format = "raster", overwrite = TRUE)

  






