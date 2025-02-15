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


# random.sp3<- generateRandomSp(my_variables, approach = "response", 
#                               realistic.sp = TRUE, 
#                               niche.breadth = 'narrow',
#                               convert.to.PA = TRUE)
# saveRDS(random.sp3, "saved_data/random_species/response_realistic_narrow.RDS")                              
#narrow makes a suitability map too small so it does not work bc it has then too many points

random.sp4<- generateRandomSp(my_variables, approach = "response", 
                              realistic.sp = TRUE, 
                              niche.breadth = 'wide',
                              convert.to.PA = TRUE)
saveRDS(random.sp4, "saved_data/random_species/response_realistic_wide.RDS")                              
 
#######################################################################
#Populations
#######################################################################
                             
#randomly select half of the popu from my random distributions
#new_random_popu<-spatSample(response_realistic_narrow$pa.raster, size = global(response_realistic_narrow$pa.raster == 1, fun = "sum", na.rm = TRUE)$sum/2 ,
                        #method= "weights", na.rm= TRUE, as.points = TRUE)
response_realistic <- readRDS("~/Desktop/UNIL/Master/Rstudio/Virtual species/saved_data/random_species/response_realistic.RDS")
response_realistic$pa.raster<-terra::unwrap(response_realistic$pa.raster)
response_realistic$suitab.raster<-terra::unwrap(response_realistic$suitab.raster)

random_total_popu<- as.points(response_realistic$pa.raster == 1, values=FALSE)
saveRDS(random_total_popu, "saved_data/random_species/random_total_popu.RDS")
sample_size<-global(response_realistic$pa.raster == 1, fun = "sum", na.rm = TRUE)$sum%/%2 
saveRDS(sample_size, "saved_data/random_species/sample_size.RDS" )
suitab_random<-response_realistic$suitab.raster

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
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 for(i in 1:repetition){
  year=year
  points(year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/stady_stady_sample/year",i,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
stady_stady_sample(6,200,random_total_popu,suitab_random )
#i put 200 for the sample size but it should be half of the number of point
#instead of 200 put sample_size
# should we do less than half should i select points from the random suitab map

#INCREASE POPU
increase_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_increase){
 #set.seed(42)
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_weights<-weights
 new_weights[indices_to_keep]<-0
 saveRDS(year, "saved_data/random_species/increase_sample/year1.RDS")
 
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
increase_sample(6,200,random_total_popu,suitab_random,45)
#DECREASE POPU
decrease_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_reduce){
 #set.seed(42)
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 inverse_weights <- 1 - weights
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_inverse_weights<-inverse_weights[indices_to_keep]
 saveRDS(year, "saved_data/random_species/decrease_sample/year1.RDS")
 
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
decrease_sample(6,200,random_total_popu,suitab_random,50)


#################################################################################
#Observations
#################################################################################

my_biais <- mixedsort(list.files(path = "biais",
                                 pattern = NULL,
                                 full.names = TRUE))
my_biais<- rast(my_biais)
head(my_biais)
names(my_biais)<-c("obs_Origin_1994", "obs_1995_2000", "obs_2001_2006", "obs_2007_2012",
                   "obs_2013_2018", "obs_2019_2024",
                   "Nobsvs_Origin_1994", "Nobsvs_1995_2000", "Nobsvs_2001_2006", "Nobsvs_2007_2012",
                   "Nobsvs_2013_2018", "Nobsvs_2019_2024",
                   "prop_pot_Origin-1994", "prop_pot_1995-2000", "prop_pot_2001-2006",  
                   "prop_pot_2007-2012", "prop_pot_2013-2018", "prop_pot_2019-2024")

#normalised the biais
biais_transformed<- my_biais
for (i in 1:length(names(my_biais))){
 biais_transformed[[i]]<-log(my_biais[[i]] + 1)
}

#count the total number of observations per raster per year
nb_obs<-data.frame(obs_Origin_1994=0, obs_1995_2000= 0, obs_2001_2006=0, obs_2007_2012=0, obs_2013_2018=0, obs_2019_2024=0)
for (i in c(1:6)){
 total_sum <- global(my_biais[[i]], fun = "sum", na.rm = TRUE)$sum
 nb_obs[,i]<-total_sum
}
#I decrease the number of sample to see if it's work otherwise it will take too long
nb_obs<-nb_obs%/%100

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
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 for(i in 1:repetition){
  year=year
  par(mfg = c(1, 1))
  points(year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/biais_stable/09_year",i,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i]], year, xy= TRUE)[, 2]
  biais_year<-sample(year, size = 0.9*nb_obs[,i], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste("saved_data/random_species/biais_stable/09_biais_year",i,".RDS", sep=""))#save each year in a file
  
 }
 return(invisible(year))
}
biais_stady_sample(6,1000,random_total_popu,suitab_random,biais_transformed)

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
 biais_year<-sample(year, size = 0.9*nb_obs[,1], prob = biais_weights, replace = TRUE)
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
  biais_year<-sample(year, size = 0.9*nb_obs[,i+1], prob = biais_weights, replace = TRUE)
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
 biais_year<-sample(year, size = 0.9*nb_obs[,1], prob = biais_weights, replace = TRUE)
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
  biais_year<-sample(year, size = 0.9*nb_obs[,i+1], prob = biais_weights, replace = TRUE)
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
 year<-sample_size[indices_to_keep,]#recreate a spatvector from the desired individuals
 for(i in 1:repetition){
  year=year
  par(mfg = c(1, 1))
  points(year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste("saved_data/random_species/biais_stable/01_year",i,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i]], year, xy= TRUE)[, 2]
  biais_year<-sample(year, size = 0.1*nb_obs[,i], prob = biais_weights, replace = TRUE)
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
 biais_year<-sample(year, size = 0.1*nb_obs[,1], prob = biais_weights, replace = TRUE)
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
  biais_year<-sample(year, size = 0.1*nb_obs[,i+1], prob = biais_weights, replace = TRUE)
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
 biais_year<-sample(year, size = 0.1*nb_obs[,1], prob = biais_weights, replace = TRUE)
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
  biais_year<-sample(year, size = 0.1*nb_obs[,i+1], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste("saved_data/random_species/biais_decrease/01_biais_year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
#biais_decrease_sample(6,popu300,sample2,prob_raster_quantile,50,biais_transformed)
