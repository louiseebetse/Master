library(virtualspecies)
library(terra)
library(geodata)
library(git2r)
library(sf)
library(gtools) # V.3.9.4
library(data.table)
library(sp)


#set the session to the right place
setwd("/Users/louise/Desktop/UNIL/Master/Rstudio/Virtual species")


#export my variables from my file
my_variables <- mixedsort(list.files(path = "variables",
                                     pattern = NULL,
                                     full.names = TRUE))
#rasterise my documents
my_variables<- rast(my_variables)
cat(crs(my_variables))
#change the names of my variables
names(my_variables) = c("ai","bio4", "bio6", "bio15","gdd3","pop2000", "pop5000", "rad", "dem", "slope", "canopy")
#compute the curvature
my_variables$dem<-terrain(my_variables$dem, v= "TPI")
#plot a raster
plot(my_variables)
#export the swiss shapefile
swiss_shape<-st_read("swiss_shapefile.gpkg")

#export my data from my file
load("SPImaster/sp/1008910")
C.kitaibelii<-clear(my.sp)
C.kitaibelii <- SpatialPointsDataFrame(coords=C.kitaibelii[,6:7], data=C.kitaibelii[ , -c(6, 7)])
C.kitaibelii<- remove.duplicates(C.kitaibelii, zero = 300)
C.kitaibelii<- as.data.frame(C.kitaibelii)
coord<- C.kitaibelii[,25:26]
#plot(C.kitaibelii)

#project the species coordinates
plot(my_variables$bio15)
points(C.kitaibelii)
variablesValues<-data.frame(extract(my_variables,coord))
head(variablesValues)

#function to clear the data
clear<- function(species){
 species$date <- as.IDate(species$date)
 species<- species[species$v_presence_status==2,]
 species<- species[species$v_introduction_status!=4,]
 species<- species[species$v_doubt_status<2,]
 species<- species[species$v_xy_radius<= 100,]
 species <- species[!is.na(species$x), ]
 species<- species[species$date>=2004,]
 species<- species[species$v_co_canton!="IT-21",]
 species<- species[species$v_co_canton!="IT-25",]
 return(species)
}



