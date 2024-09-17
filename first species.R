library(virtualspecies)
library(terra)
library(geodata)
library(git2r)
library(sf)
library(gtools) # V.3.9.4
library(data.table)
library(sp)
library(CovSel)
library(np)


#set the session to the right place
setwd("/Users/louise/Desktop/UNIL/Master/Rstudio/Virtual species")


#export my variables from my file
my_variables <- mixedsort(list.files(path = "variables",
                                     pattern = NULL,
                                     full.names = TRUE))
#rasterise my variables
my_variables<- rast(my_variables)
cat(crs(my_variables))
#change the names of my variables
names(my_variables) = c("ai","bio4", "bio6", "bio15","gdd3","pop2000", "pop5000", "rad", "dem", "slope", "canopy")
#compute the curvature
my_variables$dem<-terrain(my_variables$dem, v= "TPI")
##plot a raster
#plot(my_variables)

#export the swiss shapefile
swiss_shape<-st_read("swiss_shapefile.gpkg")
swiss_shape <- st_transform(swiss_shape, crs = 2056)

#export my species data from my file
load("SPImaster/sp/1008910")

C.kitaibelii<-clear(my.sp) #clear my species data using clear()
C.kitaibelii <- SpatialPointsDataFrame(coords=C.kitaibelii[,6:7], data=C.kitaibelii[ , -c(6, 7)]) #creating a SPDF with just the coordinates
C.kitaibelii<- remove.duplicates(C.kitaibelii, zero = 300) #remove duplicate presence in a 300 radius?
C.kitaibelii<- as.data.frame(C.kitaibelii) #come back to a data frame
coord<- C.kitaibelii[,25:26]
#plot(C.kitaibelii)


#plot(my_variables$bio15)
#points(C.kitaibelii)
#extract the variables coord at the species presence
variablesValues<-data.frame(extract(my_variables,coord))
head(variablesValues)

#put the species data in the same table
spData<- na.omit(data.frame(coord, variablesValues))
head(spData)
spData<-spData[,-3]#remove the id line
spData$presence <- rep(1, times = length(spData$x))#put a presence line

#can't use it
args<- cov.sel(T=spData$presence, Y= spData$x, X= spData[,3:13], type = "np")

glmStart<- glm(spData$presence~1, data=spData, family= "binomial")

#function to clear the data
clear<- function(species){
 species$date <- as.IDate(species$date)
 shapefile<-vect(swiss_shape$geom)
 #take out data points which are outside of the shapefile
 new_df <- data.frame(x = numeric(0), y = numeric(0))  # Create an empty data frame with same structure
 
 check_within_base <- function(i) {
  point <- vect(cbind(species$x[i], species$y[i]), type = "points")
  is_within <- relate(point, shapefile, relation = "within")
  
  if (is_within) {
   return(i)  # Return index of the row if it's within
  } else {
   return(NA)  # Return NA if not within
  }
 }
 
 # Apply function to each row index
 indices <- sapply(1:nrow(species), check_within_base)
 
 # Filter the original data frame to keep only the points that are within the shapefile
 new_df <- species[!is.na(indices), ]
 species<-new_df
 species<- species[species$v_presence_status==2,]#only presence
 species<- species[species$v_introduction_status!=4,]#
 species<- species[species$v_doubt_status<2,]#sure observation
 species<- species[species$v_xy_radius<= 100,]
 species <- species[!is.na(species$x), ]
 species<- species[species$date>=2004,]
 return(species)
}
