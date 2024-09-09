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

#export my data from my file
load("SPImaster/sp/1008910")
swiss_shape<-st_read("swiss_shapefile.gpkg")

C.kitaibelii<-clear(my.sp)
coordinates(C.kitaibelii) <- c("x", "y")
C.kitaibelii<- remove.duplicates(C.kitaibelii, zero = 300)
#plot(C.kitaibelii)

mask(swissclim,swiss_shape, inverse= FALSE)


clear<- function(species){
 species$date <- as.IDate(species$date)
 species<- species[species$v_presence_status==2,]
 species<- species[species$v_introduction_status!=4,]
 species<- species[species$v_doubt_status<2,]
 species<- species[species$v_xy_radius<= 100,]
 species <- species[!is.na(species$x), ]
 species<- species[species$date>=2004,]
 return(species)
}



