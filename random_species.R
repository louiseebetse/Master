library(virtualspecies)
library(terra)
library(gtools)
library(sf)

#set the session to the right place
setwd("/Users/louise/Desktop/UNIL/Master/Rstudio/Virtual species")

#export my variables from my file
my_variables <- mixedsort(list.files(path = "variables",
                                     pattern = NULL,
                                     full.names = TRUE))
my_variables<- rast(my_variables)
#cat(crs(my_variables))
names(my_variables) = c("ai","bio4", "bio6", "bio15","gdd3", "ph", "pop2000", "pop5000", "rad", "dem", "slope", "canopy")
my_variables$rad <- project(my_variables$rad, "EPSG:2056")
my_variables$dem<-terrain(my_variables$dem, v= "TPI")
my_variables<- subset(my_variables, c("ai","bio4", "bio6", "bio15","gdd3", "ph", "dem", "slope", "canopy")) 
my_variables<- aggregate(my_variables, fact= 4)#aggregate to have a 100m resolution

#export the swiss shapefile
swiss_shape<-st_read("swiss_shapefile.gpkg")
swiss_shape <- st_transform(swiss_shape, crs = 2056)

#generate the random species
random.sp <- generateRandomSp(my_variables, approach = "response", 
                              realistic.sp = TRUE)
#plotResponse(random.sp)
saveRDS(random.sp, "saved_data/random_species/response_realistic.RDS")                              


random.sp3<- generateRandomSp(my_variables, approach = "response", 
                              realistic.sp = TRUE, 
                              niche.breadth = 'narrow',
                              convert.to.PA = TRUE)
saveRDS(random.sp3, "saved_data/random_species/response_realistic_narrow.RDS")                              

random.sp4<- generateRandomSp(my_variables, approach = "response", 
                              realistic.sp = TRUE, 
                              niche.breadth = 'wide',
                              convert.to.PA = TRUE)
saveRDS(random.sp4, "saved_data/random_species/response_realistic_narrow.RDS")                              
 
#######################################################################
#Population and Observations
#######################################################################
                             
#randomly select half of the popu from my random distributions
#new_random_popu<-spatSample(response_realistic_narrow$pa.raster, size = global(response_realistic_narrow$pa.raster == 1, fun = "sum", na.rm = TRUE)$sum/2 ,
                        #method= "weights", na.rm= TRUE, as.points = TRUE)
random_total_popu<- as.points(response_realistic_narrow$pa.raster == 1, values=FALSE)
sample_size<-global(response_realistic_narrow$pa.raster == 1, fun = "sum", na.rm = TRUE)$sum/2 

#STABLE POPU
stady_stady_sample<- function(repetition, sample_size,sample_vector,proba_raster){
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 plot(swiss_shape$geom) #i plot the shape of switzerland
 
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition
 col <- colorRampPalette(c(start_color, end_color))(repetition)
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)# sample the individuals with weights
 year<-sample2[indices_to_keep,]#recreate a spatvector from the desired individuals
 for(i in 1:repetition){
  year=year
  points(year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/stady_stady_sample/year",i,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
stady_stady_sample(6,sample_size,new_random_popu,random.sp$pa.raster)

#INCREASE POPU
increase_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_increase){
 #set.seed(42)
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_weights<-weights
 new_weights[indices_to_keep]<-0
 saveRDS(year, "saved_data/increase_sample/year1.RDS")
 
 num_to_add <- (max_increase/(repetition-1))/100*sample_size
 
 plot(swiss_shape$geom) #i plot the shape of switzerland
 
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition 
 col <- colorRampPalette(c(start_color, end_color))(repetition) 
 points(year, col= col[1], pch= 20, cex =0.5)
 
 for(i in 1:(repetition - 1)){
  indices_to_add<-sample(nrow(sample_vector), size = num_to_add, 
                         prob = new_weights, replace = FALSE)# sample the individuals with weights
  year<-rbind(year,sample_vector[indices_to_add,])#recreate a spatvector from the desired individuals
  new_weights[indices_to_add]<-0
  points(year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/increase_sample/year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
increase_sample(6,1000,new_random_popu,random.sp$pa.raster,50)
#DECREASE POPU
decrease_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_reduce){
 #set.seed(42)
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 inverse_weights <- 1 - weights
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_inverse_weights<-inverse_weights[indices_to_keep]
 saveRDS(year, "saved_data/decrease_sample/year1.RDS")
 
 num_to_delete <- (max_reduce/(repetition-1))/100*sample_size
 
 plot(swiss_shape$geom) #i plot the shape of switzerland
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition  # Adjust as needed
 col <- colorRampPalette(c(start_color, end_color))(repetition)
 points(year, col= col[1], pch= 20, cex =0.5)
 
 for(i in 1:(repetition - 1)){
  indices_to_delete<-sample(nrow(year), size = num_to_delete, 
                            prob = new_inverse_weights, replace = FALSE)# sample the individuals with weights
  year<-year[-indices_to_delete,]#recreate a spatvector from the desired individuals
  new_inverse_weights<-new_inverse_weights[-indices_to_delete]
  points(year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/decrease_sample/year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
decrease_sample(6,1000,new_random_popu,random.sp$pa.raster,50)


#load the species 
load("SPImaster/sp/1008910")
# Convert the "date" column to Date format if it's not already
my.sp$date <- as.Date(my.sp$date)

# Define the date breaks
breaks <- as.Date(c("1994-12-31", "2000-12-31", "2006-12-31", 
                    "2012-12-31", "2018-12-31", "2024-12-31"))

# Use cut() to classify the dates into groups
categories <- cut(my.sp$date, breaks = c(as.Date("1900-01-01"), breaks), 
                  labels = c("Before 1994", "1995-2000", "2001-2006", 
                             "2007-2012", "2013-2018", "2019-2024"), 
                  right = TRUE)

#Count the number of observations in each category &
#Convert to a numeric vector
date_vector <- as.numeric(table(categories))


#STABLE OBSERVATIONS 90%
biais_stady_sample<- function(repetition, sample_size, sample_vector,proba_raster, biais_raster){
 par(mfrow = c(1, 2))
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 plot(swiss_shape$geom) #i plot the shape of switzerland
 plot(swiss_shape$geom)
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition # 6 
 col <- colorRampPalette(c(start_color, end_color))(repetition)
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)# sample the individuals with weights
 year<-sample2[indices_to_keep,]#recreate a spatvector from the desired individuals
 for(i in 1:repetition){
  year=year
  par(mfg = c(1, 1))
  points(year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/biais_stable/09_year",i,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i]], year, xy= TRUE)[, 2]
  biais_year<-sample(year, size = 0.9*date_vector[i], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste("saved_data/random_species/biais_stable/09_biais_year",i,".RDS", sep=""))#save each year in a file
  
 }
 return(invisible(year))
}
#biais_stady_sample(6,popu300,sample2,prob_raster_quantile,biais_transformed)

#INCREASE POPULATION WITH OBSERVATIONS 90%
biais_increase_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_increase,biais_raster){
 #set.seed(42)
 par(mfrow = c(1, 2))
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_weights<-weights
 new_weights[indices_to_keep]<-0
 saveRDS(year, "saved_data/random_species/biais_increase/09_year1.RDS")
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = 0.9*date_vector[1], prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year, paste("saved_data/random_species/biais_increase/09_biais_year1.RDS", sep=""))#save each year in a file
 
 
 num_to_add <- (max_increase/(repetition-1))/100*sample_size
 
 par(mfg = c(1, 1))
 plot(swiss_shape$geom) #i plot the shape of switzerland
 par(mfg = c(1, 2))
 plot(swiss_shape$geom)
 
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition 
 col <- colorRampPalette(c(start_color, end_color))(repetition)
 par(mfg = c(1, 1))
 points(year, col= col[1], pch= 20, cex =0.5)
 par(mfg = c(1, 2))
 points(unique_biais_year, col= col[1], pch= 20, cex =0.5)#make points on the map
 
 
 for(i in 1:(repetition - 1)){
  indices_to_add<-sample(nrow(sample_vector), size = num_to_add, 
                         prob = new_weights, replace = FALSE)# sample the individuals with weights
  year<-rbind(year,sample_vector[indices_to_add,])#recreate a spatvector from the desired individuals
  new_weights[indices_to_add]<-0
  par(mfg = c(1, 1))
  points(year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/biais_increase/09_year",i+1,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size+i*num_to_add, prob = biais_weights, replace = TRUE)
  biais_year<-sample(year, size = 0.9*date_vector[i+1], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste("saved_data/random_species/biais_increase/09_biais_year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
#biais_increase_sample(6,popu300,sample2,prob_raster_quantile,50,biais_transformed)

#DECREASE POPULATION WITH OBSERVATIONS 90%
biais_decrease_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_reduce,biais_raster){
 #set.seed(42)
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 inverse_weights <- 1 - weights
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_inverse_weights<-inverse_weights[indices_to_keep]
 saveRDS(year, "saved_data/random_species/biais_decrease/09_year1.RDS")
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = 0.9*date_vector[1], prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year, paste("saved_data/random_species/biais_decrease/09_biais_year1.RDS", sep=""))#save each year in a file
 
 
 num_to_delete <- (max_reduce/(repetition-1))/100*sample_size
 
 par(mfg = c(1, 1))
 plot(swiss_shape$geom) #i plot the shape of switzerland
 par(mfg = c(1, 2))
 plot(swiss_shape$geom)
 
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition 
 col <- colorRampPalette(c(start_color, end_color))(repetition)
 par(mfg = c(1, 1))
 points(year, col= col[1], pch= 20, cex =0.5)
 par(mfg = c(1, 2))
 points(unique_biais_year, col= col[1], pch= 20, cex =0.5)#make points on the map
 
 for(i in 1:(repetition - 1)){
  indices_to_delete<-sample(nrow(year), size = num_to_delete, 
                            prob = new_inverse_weights, replace = FALSE)# sample the individuals with weights
  year<-year[-indices_to_delete,]#recreate a spatvector from the desired individuals
  new_inverse_weights<-new_inverse_weights[-indices_to_delete]
  par(mfg = c(1, 1))
  points(year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/biais_decrease/09_year",i+1,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size-i*num_to_delete, prob = biais_weights, replace = TRUE)
  biais_year<-sample(year, size = 0.9*date_vector[i+1], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste("saved_data/random_species/biais_decrease/09_biais_year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
#biais_decrease_sample(6,popu300,sample2,prob_raster_quantile,50,biais_transformed)


#STABLE OBSERVATIONS 10%
biais_stady_sample<- function(repetition, sample_size, sample_vector,proba_raster, biais_raster){
 par(mfrow = c(1, 2))
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 plot(swiss_shape$geom) #i plot the shape of switzerland
 plot(swiss_shape$geom)
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition # 6 
 col <- colorRampPalette(c(start_color, end_color))(repetition)
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)# sample the individuals with weights
 year<-sample2[indices_to_keep,]#recreate a spatvector from the desired individuals
 for(i in 1:repetition){
  year=year
  par(mfg = c(1, 1))
  points(year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/biais_stable/01_year",i,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i]], year, xy= TRUE)[, 2]
  biais_year<-sample(year, size = 0.1*date_vector[i], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste("saved_data/random_species/biais_stable/01_biais_year",i,".RDS", sep=""))#save each year in a file
  
 }
 return(invisible(year))
}
#biais_stady_sample(6,popu300,sample2,prob_raster_quantile,biais_transformed)

#INCREASE POPULATION WITH OBSERVATIONS 10%
biais_increase_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_increase,biais_raster){
 #set.seed(42)
 par(mfrow = c(1, 2))
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_weights<-weights
 new_weights[indices_to_keep]<-0
 saveRDS(year, "saved_data/random_species/biais_increase/01_year1.RDS")
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = 0.1*date_vector[1], prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year, paste("saved_data/random_species/biais_increase/01_biais_year1.RDS", sep=""))#save each year in a file
 
 
 num_to_add <- (max_increase/(repetition-1))/100*sample_size
 
 par(mfg = c(1, 1))
 plot(swiss_shape$geom) #i plot the shape of switzerland
 par(mfg = c(1, 2))
 plot(swiss_shape$geom)
 
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition 
 col <- colorRampPalette(c(start_color, end_color))(repetition)
 par(mfg = c(1, 1))
 points(year, col= col[1], pch= 20, cex =0.5)
 par(mfg = c(1, 2))
 points(unique_biais_year, col= col[1], pch= 20, cex =0.5)#make points on the map
 
 
 for(i in 1:(repetition - 1)){
  indices_to_add<-sample(nrow(sample_vector), size = num_to_add, 
                         prob = new_weights, replace = FALSE)# sample the individuals with weights
  year<-rbind(year,sample_vector[indices_to_add,])#recreate a spatvector from the desired individuals
  new_weights[indices_to_add]<-0
  par(mfg = c(1, 1))
  points(year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/biais_increase/01_year",i+1,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size+i*num_to_add, prob = biais_weights, replace = TRUE)
  biais_year<-sample(year, size = 0.1*date_vector[i+1], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste("saved_data/random_species/biais_increase/01_biais_year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
#biais_increase_sample(6,popu300,sample2,prob_raster_quantile,50,biais_transformed)

#DECREASE POPULATION WITH OBSERVATIONS 10%
biais_decrease_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_reduce,biais_raster){
 #set.seed(42)
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 inverse_weights <- 1 - weights
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_inverse_weights<-inverse_weights[indices_to_keep]
 saveRDS(year, "saved_data/random_species/biais_decrease/01_year1.RDS")
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = 0.1*date_vector[1], prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year, paste("saved_data/random_species/biais_decrease/01_biais_year1.RDS", sep=""))#save each year in a file
 
 
 num_to_delete <- (max_reduce/(repetition-1))/100*sample_size
 
 par(mfg = c(1, 1))
 plot(swiss_shape$geom) #i plot the shape of switzerland
 par(mfg = c(1, 2))
 plot(swiss_shape$geom)
 
 start_color <- "#1c74c1" # First color
 end_color <- "#fff064"   # Last color
 repetition <- repetition 
 col <- colorRampPalette(c(start_color, end_color))(repetition)
 par(mfg = c(1, 1))
 points(year, col= col[1], pch= 20, cex =0.5)
 par(mfg = c(1, 2))
 points(unique_biais_year, col= col[1], pch= 20, cex =0.5)#make points on the map
 
 for(i in 1:(repetition - 1)){
  indices_to_delete<-sample(nrow(year), size = num_to_delete, 
                            prob = new_inverse_weights, replace = FALSE)# sample the individuals with weights
  year<-year[-indices_to_delete,]#recreate a spatvector from the desired individuals
  new_inverse_weights<-new_inverse_weights[-indices_to_delete]
  par(mfg = c(1, 1))
  points(year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/biais_decrease/01_year",i+1,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size-i*num_to_delete, prob = biais_weights, replace = TRUE)
  biais_year<-sample(year, size = 0.1*date_vector[i+1], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste("saved_data/random_species/biais_decrease/01_biais_year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
#biais_decrease_sample(6,popu300,sample2,prob_raster_quantile,50,biais_transformed)
