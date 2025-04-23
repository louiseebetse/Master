library(virtualspecies)
library(terra)
library(gtools)
library(sf)
library(geodata)
library(git2r)
library(data.table)
library(sp)
library(biomod2)
library(MASS)
library(reshape2)
library(tidyr)
library(readxl)
library(httr)
library(mapview)
library(plot.matrix)
library(magick)
library(dplyr)
library(raster)


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


biais_stady_sample<- function(repetition, sample_size, sample_vector,proba_raster, biais_raster,path, percentage_obs){
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
  saveRDS(year, paste(path, "/year",i,".RDS", sep=""))#save each year in a file
  
  
  biais_weights<- terra::extract(biais_raster[[i]], year, xy= TRUE)[, 2]

  
  if (as.integer(percentage_obs[2,i]*sample_size)==0){ #add a condition if there are 0 observations this year
   unique_biais_year<-vect() #there will be 0 observation this year
  }
  else{
   biais_year<-sample(year, size = as.integer(percentage_obs[2,i]*sample_size), prob = biais_weights, replace = TRUE)
   unique_biais_year <- terra::unique(biais_year)
  }
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste(path, "/biais_year",i,".RDS", sep=""))#save each year in a file
  
  
 }
 return(invisible(year))
}
biais_increase_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_increase,biais_raster,path,percentage_obs){
 # #set.seed(42)
 # repetition= 14
 # sample_size= 100
 # sample_vector= sample2
 # proba_raster= prob_raster_quantile
 # max_increase=50
 # biais_raster= biais_transformed
 # path= path_1
 # obs_by_year= date_vector
 # 
 #  
 
 par(mfrow = c(1, 2))
 
 start_weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 starting_points<- sample(nrow(sample_vector), size = 3, #create the first 3 points from where the popu will spread
                          prob = start_weights, replace = FALSE)
 starting_points<-sample_vector[starting_points,]
 radius_rast<-ifel(is.na(proba_raster), NA, 0)# create empty raster with same resolution #fill it with 0 except where there are NAs
 dist_raster <- distance(radius_rast, starting_points)  # Compute distance raster
 d50 <- 8000  # Distance where probability is 0.5
 k <- -0.001  # Controls steepness of the decay
 
 dist_raster <-(1-( 1 / (1 + exp(k * (values(dist_raster) - d50))))) #0.5 in a 5km radius of the point then decay
 dist_raster <- (dist_raster - min(dist_raster)) / (max(dist_raster) - min(dist_raster)) * 0.3
 values(radius_rast) <- dist_raster  # Assign computed probabilities
 radius_rast <- mask(radius_rast, proba_raster) #put the NA to keep the shape
 par(mfg = c(1, 1))
 #plot(radius_rast)
 
 weights <- terra::extract(radius_rast, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 #new_weights<-weights
 #new_weights[indices_to_keep]<-0
 saveRDS(year, paste(path, "/year1.RDS", sep="")) #save each year in a file
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = as.integer(percentage_obs[2,1]*sample_size), prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year, paste(path, "/biais_year1.RDS", sep=""))
 
 num_to_add <- (max_increase/(repetition-1))/100*sample_size
 if (num_to_add <1) {num_to_add<-1} #add at least one popu
 
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
  
  radius_rast<-ifel(is.na(proba_raster), NA, 0)# create empty raster with same resolution #fill it with 0 except where there are NAs
  dist_raster <- distance(radius_rast, starting_points)  # Compute distance raster
  d50 <- d50+12000  # Distance where probability is 0.5
  k <- -0.0001  # Controls steepness of the decay
  
  dist_raster <-(1-( 1 / (1 + exp(k * (values(dist_raster) - d50))))) #0.5 in a 5km radius of the point then decay
  dist_raster <- (dist_raster - min(dist_raster)) / (max(dist_raster) - min(dist_raster)) * 0.3
  values(radius_rast) <- dist_raster  # Assign computed probabilities
  
  # Assign computed probabilities
  radius_rast <- mask(radius_rast, proba_raster) #put the NA to keep the shape
  par(mfg = c(1, 1))
  #plot(radius_rast)
  
  new_weights <- terra::extract(radius_rast, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
  indices_to_add<-sample(nrow(sample_vector), size = num_to_add, 
                         prob = new_weights, replace = FALSE)# sample the individuals with weights
  year<-rbind(year,sample_vector[indices_to_add,])#recreate a spatvector from the desired individuals
  #new_weights[indices_to_add]<-0
  par(mfg = c(1, 1))
  points(year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(year, paste(path, "/year",i+1,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size+i*num_to_add, prob = biais_weights, replace = TRUE)
  if (as.integer(percentage_obs[2,i+1]*sample_size)==0){
   unique_biais_year<-NA
  }
  else{
   biais_year<-sample(year, size = as.integer(percentage_obs[2,i+1]*sample_size), prob = biais_weights, replace = TRUE)
   unique_biais_year <- terra::unique(biais_year)
  }
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste(path, "/biais_year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
biais_decrease_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_reduce,biais_raster,path,percentage_obs){
 #set.seed(42)
 par(mfrow = c(1, 2))
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 inverse_weights <- 1 - weights
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_inverse_weights<-inverse_weights[indices_to_keep]
 saveRDS(year, paste(path, "/year1.RDS", sep=""))
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = as.integer(percentage_obs[2,1]*sample_size), prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year, paste(path, "/biais_year1.RDS", sep=""))
 num_to_delete <- (max_reduce/(repetition-1))/100*sample_size
 if (num_to_delete <1) {num_to_delete<-1} #delete at least one population
 
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
  saveRDS(year, paste(path, "/year",i+1,".RDS", sep=""))#save each year in a file
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size-i*num_to_delete, prob = biais_weights, replace = TRUE)
  if (as.integer(percentage_obs[2,i+1]*sample_size)==0){
   unique_biais_year<-vect()
  }
  else{
   biais_year<-sample(year, size = as.integer(percentage_obs[2,i+1]*sample_size), prob = biais_weights, replace = TRUE)
   unique_biais_year <- terra::unique(biais_year)
  }
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste(path, "/biais_year",i+1,".RDS", sep=""))#save each year in a file
  
 }
 return(invisible(year)) 
}


#################################################################################
#Populations size
#################################################################################
data_list <- list()
# Loop through each file
for (file in list.files(path = "SPImaster/sp", full.names = TRUE)) {
 # Read the file
 load(file)
 data<-my.sp
 file_name <- tools::file_path_sans_ext(basename(file))
 # Store the data in the list with the file name as the key
 data_list[[file_name]] <- data
}

#get the number of observations by plant after cleaning
my_size_popu<- function(list_presences){
 #list_distribution<- list()
 name<- names(list_presences)
 #store the size of the popus after desagreging
 size_popus<- c(1:length(list_presences))
 for (sp in seq_along(list_presences)) {
  par(mfrow = c(1, 1))
  ## CLEAR THE DATA
  
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
  
  species<-clear(list_presences[[sp]]) #clear my species data using clear()
  species <- SpatialPointsDataFrame(coords=species[,c("x", "y")], data=species[ , -c(6, 7)]) #creating a SPDF with just the coordinates
  species<- remove.duplicates(species, zero = 500) #remove duplicate presence in a 300 radius?
  species<- as.data.frame(species) #come back to a data frame
  coord<- species[,c("x", "y")]
  size_popus[sp]<-as.numeric(length(coord$x))
 }
 return(size_popus)
}
popu_size<-my_size_popu(data_list)

mean_popu_size<- mean(popu_size) #get the mean of the popu_size
sd_popu_size <- sd(popu_size) #get the sd of the popu_size

#randomly draw a value
#popu_size<-as.integer(sqrt(rnorm(1, mean=mean_popu_size, sd=sd_popu_size)^2))
# test <- -1  # initialize to a negative value
# while (test <= 0) { test <- as.integer(rnorm(1, mean = mean_popu_size, sd = sd_popu_size))}
#################################################################################
#Observations
#################################################################################

my_biais <- mixedsort(list.files(path = "biais",
                                 pattern = NULL,
                                 full.names = TRUE))
my_biais<- rast(my_biais)
head(my_biais)
names(my_biais)<-c("obs_Origin_1946", "obs_1947_1952", "obs_1953_1958", "obs_1959_1964", "obs_1965_1970", 
                   "obs_1971_1976","obs_1977_1982", "obs_1983_1988","obs_1989_1994", "obs_1995_2000", 
                   "obs_2001_2006", "obs_2007_2012","obs_2013_2018", "obs_2019_2024",
                   "Nobsv_Origin_1946", "Nobsv_1947_1952", "Nobsv_1953_1958", "Nobsv_1959_1964", "Nobsv_1965_1970", 
                   "Nobsv_1971_1976","Nobsv_1977_1982", "Nobsv_1983_1988","Nobsv_1989_1994", "Nobsv_1995_2000", 
                   "Nobsv_2001_2006", "Nobsv_2007_2012","Nobsv_2013_2018", "Nobsv_2019_2024",
                   "prop_pot_1977_1982", "prop_pot_1947_1952","prop_pot_1953_1958", "prop_pot_1959_1964", 
                   "prop_pot_1965_1970", "prop_pot_1971_1976")

#normalised the biais
biais_transformed<- my_biais
for (i in 1:length(names(my_biais))){
 biais_transformed[[i]]<-log(my_biais[[i]] + 1)
}

#count the total number of observations per raster per year
nb_obs<-data.frame(obs_Origin_1946= c(0,0), obs_1947_1952= c(0,0), obs_1953_1958= c(0,0),
                   obs_1959_1964= c(0,0), obs_1965_1970= c(0,0), 
                   obs_1971_1976= c(0,0), obs_1977_1982 = c(0,0), obs_1983_1988= c(0,0),
                   obs_1989_1994=c(0,0), obs_1995_2000= c(0,0), obs_2001_2006=c(0,0),
                   obs_2007_2012=c(0,0), obs_2013_2018=c(0,0), obs_2019_2024=c(0,0))

#add the total number of obs per year in the first column of the df
for (i in c(1:14)){
 total_sum <- global(my_biais[[i]], fun = "sum", na.rm = TRUE)$sum
 nb_obs[1,i]<-total_sum
}
total_obs<-sum(nb_obs)
#put the percentage of obs on the total number of obs all year counted
for (i in c(1:14)){
 percentage <- nb_obs[1,i]/total_obs
 nb_obs[2,i]<-percentage
}

#################################################################################
#random distributions
#################################################################################
nb_random_species <- 10
for (j in 1:nb_random_species){
 dir.create(file.path(paste("saved_data/random_species/random_sp_", j,sep="")), showWarnings = FALSE)
 random.sp4<- generateRandomSp(my_variables, approach = "response", 
                              realistic.sp = TRUE, 
                              niche.breadth = 'wide',
                              convert.to.PA = TRUE)
print(paste("random species nb",j,"has been generated"))
#######################################################################
#Populations
#######################################################################
#I remove the points that were in a probability lower than 20% to have less populations
#random.sp4$pa.raster[random.sp4$suitab.raster<0.2]<-0
suitab_random<-random.sp4$suitab.raster
suitab_random[suitab_random<0.05]<-0

#randomly select a number for the popu
nb_popu <- -1  # initialize to a negative value
while (nb_popu <= 0) { nb_popu <- as.integer(rnorm(1, mean = mean_popu_size, sd = sd_popu_size))}#loop until finds a positive value

#sample the popus based on the random number of popu *3 and weigthed by suitability raster
sample2<- spatSample(suitab_random,size= nb_popu*3, method= "weights", na.rm= TRUE, as.points = TRUE)
#size_popus= nb_popu


#random_total_popu<- as.points(nb_popu, values=FALSE)
#saveRDS(random_total_popu, "saved_data/random_species/random_total_popu.RDS")
#sample_size<-global(random.sp4$pa.raster == 1, fun = "sum", na.rm = TRUE)$sum%/%3
if (length(sample2)<=nb_popu){nb_popu=length(sample2)%/%3}
else{nb_popu<- nb_popu}

#saveRDS(sample_size, "saved_data/random_species/sample_size.RDS" )


print(paste("the population nb",j,"has been generated"))

#I decrease the number of sample to see if it's work otherwise it will take too long
#nb_obs<-nb_obs%/%100


#STABLE OBSERVATIONS
path_1= paste("saved_data/random_species/random_sp_",j,"/biais_stable", sep="")
dir.create(file.path(path_1), showWarnings = FALSE)
biais_stady_sample(repetition=14,nb_popu,sample2,suitab_random,biais_transformed, path_1, nb_obs)


#INCREASE POPULATION WITH OBSERVATIONS 50%
path_1= paste("saved_data/random_species/random_sp_",j,"/biais_increase", sep="")
dir.create(file.path(path_1), showWarnings = FALSE)
biais_increase_sample(repetition=14,nb_popu,sample2,suitab_random,80,biais_transformed, path_1, nb_obs)

#DECREASE POPULATION WITH OBSERVATIONS 50%
path_1= paste("saved_data/random_species/random_sp_",j,"/biais_decrease", sep="")
dir.create(file.path(path_1), showWarnings = FALSE)
biais_decrease_sample(repetition=14,nb_popu,sample2,suitab_random,80,biais_transformed,path_1, nb_obs)

print(paste("stable + increase + decrease of random sp",j,"has been generated"))

saveRDS(suitab_random, paste("saved_data/random_species/random_sp_", j,"/biais_stable/suitab_random.RDS",sep=""))
saveRDS(suitab_random, paste("saved_data/random_species/random_sp_", j,"/biais_increase/suitab_random.RDS",sep=""))
saveRDS(suitab_random, paste("saved_data/random_species/random_sp_", j,"/biais_decrease/suitab_random.RDS",sep=""))
}


#################################################################################
#5km distribtion
#################################################################################

as.square<- function(popu,obs){
 popu<- terra::unwrap(popu)
 obs<- terra::unwrap(obs)
 if (nrow(obs)==0){
  counts_raster_biais<-ifel(is.na(my_biais$obs_1947_1952), NA, 0)# create empty raster with same resolution #fill it with 0 except where there are NAs
 }
 else{
  biais_cell_ids <- cellFromXY(my_biais$obs_1947_1952, crds(obs))
  counts_biais <- table(biais_cell_ids)
  counts_raster_biais <- my_biais$obs_1947_1952
  values(counts_raster_biais) <- ifelse(is.na(values(counts_raster_biais)), NA, 0) # Initialize all cells with 0
  values(counts_raster_biais)[as.numeric(names(counts_biais))] <- counts_biais
 }
 #### put back in a 5km raster 
 # Extract cell IDs for each observation
 cell_ids <- cellFromXY(my_biais$obs_1947_1952, crds(popu)) 
 #obs_Origin n'a pas besoin d'etre changer d'une iteration à l'autre
 
 # Count the number of observations per cell
 counts_cell <- table(cell_ids)
 #ICI on sait deja la taille
 
 # Create an empty raster with the same structure as your original raster
 counts_raster <- my_biais$obs_1947_1952
 values(counts_raster) <-ifelse(is.na(values(counts_raster)), NA, 0)  # Initialize all cells with 0
 
 # Assign counts to the appropriate cells
 values(counts_raster)[as.numeric(names(counts_cell))] <- counts_cell
 
 # Plot the raster of counts
 #plot(counts_raster, main = "Population Counts per 5km Pixel")
 #plot(counts_raster_biais, main = "Observation Counts per 5km Pixel")
 
 #test[[1]]= counts_raster
 #test[[2]]= counts_raster_biais
 return( c(counts_raster, counts_raster_biais))
 
}

species_files<-list.files(path = "saved_data/random_species", pattern = NULL,full.names = TRUE)
for (sp in 1:(length(species_files))){
 dir.create(file.path(paste(species_files[sp],"/sq_biais_stable",sep="")), showWarnings = FALSE)
 dir.create(file.path(paste(species_files[sp],"/sq_biais_increase",sep="")), showWarnings = FALSE)
 dir.create(file.path(paste(species_files[sp],"/sq_biais_decrease",sep="")), showWarnings = FALSE)
 repet <- length(list.files(path = paste(species_files[sp],"/biais_increase/", sep="")))%/%2
 for (yr in 1:repet){
  #stable
  popu_stable<- readRDS(paste(species_files[sp],"/biais_stable/year", yr, ".RDS", sep=""))
  obs_stable<-readRDS(paste(species_files[sp],"/biais_stable/biais_year", yr, ".RDS", sep=""))
  squares<- as.square(popu_stable,obs_stable)
  saveRDS(squares[[1]], paste(species_files[sp],"/sq_biais_stable/year",yr,".RDS", sep=""))
  saveRDS(squares[[2]], paste(species_files[sp],"/sq_biais_stable/biais_year",yr,".RDS", sep=""))
  
  #increase
  popu_increase<- readRDS(paste(species_files[sp],"/biais_increase/year", yr, ".RDS", sep=""))
  obs_increase<-readRDS(paste(species_files[sp],"/biais_increase/biais_year", yr, ".RDS", sep=""))
  squares<- as.square(popu_increase,obs_increase)
  #print(paste0(yr, " ",1.5))
  saveRDS(squares[[1]], paste(species_files[sp],"/sq_biais_increase/year",yr,".RDS", sep=""))
  saveRDS(squares[[2]], paste(species_files[sp],"/sq_biais_increase/biais_year",yr,".RDS", sep=""))
  
  #decrease
  #print(paste0(yr, 3))
  popu_decrease<- readRDS(paste(species_files[sp],"/biais_decrease/year", yr, ".RDS", sep=""))
  obs_decrease<-readRDS(paste(species_files[sp],"/biais_decrease/biais_year", yr, ".RDS", sep=""))
  squares<- as.square(popu_decrease,obs_decrease)
  saveRDS(squares[[1]], paste(species_files[sp],"/sq_biais_decrease/year",yr,".RDS", sep=""))
  saveRDS(squares[[2]], paste(species_files[sp],"/sq_biais_decrease/biais_year",yr,".RDS", sep=""))
  
 }
 print(paste0(species_files[sp], " is aggregated"))
} 

#################################################################################
#table trend
#################################################################################


#### table of 5km square for popu and observations
trends<- c("biais_stable", "biais_increase", "biais_decrease" )
#acess the data saved and create table of popu +obs
#= "saved_data/random_species/random_sp_2/biais_increase"
table_maker<- function(raster5x5 = my_biais$obs_1947_1952, file.path ){
 all_the_files<-list.files(path = file.path,
                           pattern = NULL,
                           full.names = TRUE)
 # Extract the number after "year" using a regular expression
 year_numbers <- as.numeric(gsub(".*year([0-9]+)\\.RDS", "\\1", all_the_files))
 
 # Reorder the files by those extracted numbers
 all_the_files <- all_the_files[order(year_numbers)]
 
 data <- matrix(NA, nrow = 4, ncol = length(all_the_files)/2) #create the matrix
 colnames(data) <- paste0("TimeStep", 1:14) #put the name of the columns
 rownames(data) <- c("Biais_sq", "Population_sq", "Biais", "Population") #put names of the lines
 data <- as.data.frame(data)
 for (i in 1:(length(all_the_files)/2)){ #loop for the number of year
  biais_year<-readRDS(all_the_files[i*2-1])#load the observations
  year<-readRDS(all_the_files[i*2])#load the populations
  
  biais_cell_ids <- cellFromXY(raster5x5, crds(biais_year)) #aggregate to 5km
  cell_ids <- cellFromXY(raster5x5, crds(year)) #aggregate to 5km
  
  # Count the number of observations per cell
  counts_cell <- length(table(cell_ids)) #count the number of cells at 5km
  counts_biais <- length(table(biais_cell_ids)) #count the number of cells at 5km
  
  data[1, i] <- counts_biais #nb of square of observations
  data[2, i] <- counts_cell #nb of square of populations
  data[3, i] <- length(biais_year)# nb of observations
  data[4, i] <- length(year)# nb of populations
 }
 return(data)
}

# Create a plot from the dataframe created in decrease table
#add Observations
trend_maker<- function(table,species,trend){
 date<-c(1946, 1952, 1958, 1964, 1970, 1976, 1982, 1988, 1994,2000,2006,2012,2018,2024)
 id <- sub(".*sp_", "", species)
 plot(date, decrease_table[3,], type = "o", col = "#359B73", 
      xlab = "Time", ylab = "Nb of points", 
      xlim = c(1946,2024), ylim = c(0,1.8*max(decrease_table)), 
      main = paste(trend, "of species", id,  "over Time"), lty= "44",
      xaxt = "n")
 axis(1, at = date, las=2)
 # Add the second line (Population)
 lines(date, decrease_table[4,], type = "o", col = "#359B73")
 
 # Add the sq Observations
 lines(date, decrease_table[1,], type = "o", col = "#3DB7E9",lty= "44")
 
 # Add the sq Populations
 lines(date, decrease_table[2,], type = "o", col = "#3DB7E9")
 
 # Add a legend
 legend("topright", legend = c("100m Popu", "100m Obs", "5000m Popu", "5000m Obs"), 
        col = c("#359B73", "#359B73", "#3DB7E9", "#3DB7E9"), lty= c( "solid", "44", "solid","44"), pch = 1, ncol =2)
 
}

#loop through the species files to create the table and the graph to all species with their trend
#sp=2
for (sp in 1:(length(species_files))){ #loop through random species
 dir.create(file.path(paste(species_files[sp],"/table_trend", sep="")), showWarnings = FALSE)
 id <- sub(".*sp_", "", species_files[sp])
 pdf(paste(species_files[sp],"/table_trend/graph_trend_", id,".pdf",sep=""), width = 6.4, height = 7.26)
 par(mfrow = c(3, 1))
 
 for (tr in 1:(length(trends))){ #loop through the three trend
  decrease_table<-table_maker(file.path = paste(species_files[sp],"/", trends[tr],sep="")) #apply the function to create the table
  par(mfg = c(tr, 1))
  trend_maker(decrease_table,species_files[sp],trends[tr]) # create the graph for for the trend now
  saveRDS(decrease_table, paste(species_files[sp],"/table_trend/",trends[tr],".RDS",sep=""))
 }
 dev.off() 
}


#################################################################################
#Occupancy parameters
#################################################################################

# parameters
res <- 5000
pct_ws_in_5x5_threshold<-30 # should be replaced by the area * the rarity * the suitability

# load main data
### Welten and Sutter to 5x5 km square
ws_pct <-
 read.csv("data/ws_5x5_pct.txt", sep = ";") # percent of WS atlas on square from HES (check documentation/data_info.rtf)

### Species strategy and informations (check documentation/data_info.rtf)
strategy <-
 read_excel("data/SPI_TaxaInfos_v2_20220810.xlsx") # taxa info compilation
first_year_mention <- read.csv("data/year_first_mention.csv", sep = ";")
ch <- st_read(
 'data/chBoundary/SHAPEFILE_LV95_LN02/swissBOUNDARIES3D_1_4_TLM_LANDESGEBIET.shp'
) %>%
 filter(ICC %in% c('CH', "LI")) %>%
 # st_geometry %>%
 st_transform(2056) #swiss mask
ch_ext <- st_bbox(ch)

ch_sq <- st_make_grid(ch, res, offset = c(floor(ch_ext / res)[c(1, 2)] *
                                           res))[ch] #create grid aligned on 5km multpiles
sq_id <- st_coordinates(st_centroid(ch_sq)) - (res / 2) # get the coordinates of the left down corner
sq_id[, 1] <- sq_id[, 1] - 2000000 # remove the ch1903+ corrections
sq_id[, 2] <- sq_id[, 2] - 1000000  # remove the ch1903+ corrections
sq_id <- sq_id / (10 ^ (floor(log10(res)))) # round to the multiple  of the closest power of ten
sq_id <- sq_id[, 1] * (10 ^ (floor(log10(res)))) + sq_id [, 2] # merge coordinates to get the spatial id

ch_sq <- st_as_sf(data.frame(sq_id = sq_id, geometry = ch_sq))

#load the file with the biais raster
my_biais <- mixedsort(list.files(path = "biais",
                                 pattern = NULL,
                                 full.names = TRUE))
my_biais<- rast(my_biais)
my_biais<-terra::unwrap(my_biais)
#head(my_biais)
names(my_biais)<-c("obs_Origin_1946", "obs_1947_1952", "obs_1953_1958", "obs_1959_1964", "obs_1965_1970", 
                   "obs_1971_1976","obs_1977_1982", "obs_1983_1988","obs_1989_1994", "obs_1995_2000", 
                   "obs_2001_2006", "obs_2007_2012","obs_2013_2018", "obs_2019_2024",
                   "Nobsv_Origin_1946", "Nobsv_1947_1952", "Nobsv_1953_1958", "Nobsv_1959_1964", "Nobsv_1965_1970", 
                   "Nobsv_1971_1976","Nobsv_1977_1982", "Nobsv_1983_1988","Nobsv_1989_1994", "Nobsv_1995_2000", 
                   "Nobsv_2001_2006", "Nobsv_2007_2012","Nobsv_2013_2018", "Nobsv_2019_2024",
                   "prop_pot_1977_1982", "prop_pot_1947_1952","prop_pot_1953_1958", "prop_pot_1959_1964", 
                   "prop_pot_1965_1970", "prop_pot_1971_1976")

#apply the log scale to my biais raster
biais_transformed<- my_biais
for (i in 1:length(names(my_biais))){
 biais_transformed[[i]]<-log(my_biais[[i]] + 1)
}

#################################################################################
#Occupancy matrix
#################################################################################

list_species<-list.files(path = paste0("saved_data/random_species/"), full.names = FALSE)
pt_list<-list_species[4:17]
#pt_list<-pt_list[2]
trends<- c("biais_stable", "biais_increase", "biais_decrease" )
timescale<-c(1946, 1952, 1958, 1964, 1970, 1976, 1982, 1988, 1994,2000,2006,2012,2018,2024)
#path= "saved_data/species/sp_1008910/biais_increase/biais_year4.RDS" 

to_matrix<-function(path){
 id<- as.numeric(sub(".*random_sp_([0-9]+)/.*", "\\1", path))
 #fix the name
 year<- as.numeric(gsub(".*year([0-9]+)\\.RDS", "\\1", path)) #fix the year
 #print(year)
 points<-readRDS(path)  # my observations
 polygons <- ch_sq # Load polygon data
 
 points<-sf::st_as_sf(points)
 #----------------------------------------------------------------------------------
 #suitability 100m in 5km resolution
 
 suitab_map <- prob_raster_quantile
 suitab_df<- data.frame(sq_id=ch_sq$sq_id)
 
 spatvector <- vect(ch_sq)  # Convert 5km poly to spatvector
 mean_extract <- terra::extract(suitab_map, spatvector, fun = mean, na.rm = TRUE)
 suitab_df$suitability <- mean_extract[, 2]
 #----------------------------------------------------------------------------------
 #sampling effort
 nb_obs<- biais_transformed[[1:14]]
 obs_extract<- terra::extract(nb_obs[[year]],spatvector)
 suitab_df$sampling_eff <- obs_extract[,2]
 #----------------------------------------------------------------------------------
 # Perform spatial intersection: Find which points fall within which polygon
 
 #if more than 1 obs
 if (nrow(points)!=0){
  intersections <- sf::st_join(points, polygons)
  
  # Count the number of points per polygon
  point_counts <- raster::aggregate(intersections, by = list(intersections$sq_id), FUN = length)
  
  # Convert the result into a data frame
  result_df <- data.frame(
   sq_id = point_counts$Group.1,
   num_points = point_counts$sq_id  
  )
  # Merge nb of observation and suitability
  my_df <- merge(suitab_df, result_df, by = "sq_id", all.x = TRUE)
 }
 
 #if 0 obs
 else{
  my_df<-suitab_df
  my_df[, 'num_points']<-0
 }
 #----------------------------------------------------------------------------------
 
 
 # Replace NA values in 'value2' with 0
 my_df$num_points[is.na(my_df$num_points)] <- 0
 my_df$suitability[is.na(my_df$suitability)] <- 0
 my_df$year <- timescale[year]
 my_df$species <- id
 
 return(my_df)
}
occupancy_function2 <- function(sp.info, strategy_matrix,sp_matrix, effort_matrix,suit_matrix){
 # Generate sigmoid curves for colonization and persistence ---------------
 # For colonization
 t_seq <- seq(0, sp.info$colonisation_time - 1, 1)
 k_val <- 8/ length(t_seq)   # Adjust steepness (higher k = sharper transition)
 t0 <- max(t_seq) / 2         # Midpoint of the transition
 colo_coef <- 1 / (1 + exp(-k_val * (t_seq - t0)))
 colo_coef<- append(colo_coef, 1) #add 1 as the observation at the end of the vector
 #if the probability are smaller than the nb of time step i add 0 to it 
 if (length(colo_coef)<length(timescale)){colo_coef <- c(rep(0, length(timescale) - length(colo_coef)), colo_coef)
 } 
 
 
 
 # For persistence
 t_seq <- seq(0, sp.info$persistence_time - 1, 1)
 k_val <- 8 / length(t_seq)
 t0 <- max(t_seq) / 2
 persi_coef <- 1 / (1 + exp(-k_val * (t_seq - t0)))
 persi_coef<- append(persi_coef, 1)
 if (length(persi_coef)<length(timescale)){persi_coef <- c(rep(0, length(timescale) - length(persi_coef)), persi_coef)} 
 
 # Correction of occupancy probability based on colonization hypothesis
 col_coef <- sp_matrix * 0  # Initialize with 0
 per_coef <- sp_matrix * 0
 sam_coef <- sp_matrix * 0
 neo_coef <- sp_matrix * 0  # Add neophyte term if necessary
 
 # colonization matrix
 for (i in 1:nrow (sp_matrix)){
  obsi <- which(sp_matrix[i,]== 1)
  tj <- 1
  for (j in 1:length(obsi)){
   nyearj <- min(obsi[j]-tj,length(colo_coef))
   col_coefj <- colo_coef[(length(colo_coef) - nyearj):length(colo_coef)]
   col_coef[i,tj:obsi[j]] <- col_coefj
   tj <- obsi[j]+1
  }
 }
 
 
 # # Persistence matrix
 for (i in 1:nrow(sp_matrix)) {
  obsi <- which(sp_matrix[i, ] == 1)
  for (j in seq_along(obsi)) {
   nyearj <- min(ncol(sp_matrix) - obsi[j] + 1, length(persi_coef))
   if(nyearj > 0){
    per_coefj <- rev(persi_coef)[1:nyearj]
    per_coef[i, obsi[j]:(obsi[j] + nyearj - 1)] <- per_coefj
   }
  }
 }
 
 # # Colonization matrix
 # for (i in 1:nrow(sp_matrix)) {
 #  obsi <- which(sp_matrix[i, ] == 1)
 #  tj <- 1
 #  for (j in seq_along(obsi)) {
 #   nyearj <- min(obsi[j] - tj, length(colo_coef))
 #   if(nyearj > 0){
 #    col_coefj <- colo_coef[(length(colo_coef) - nyearj):length(colo_coef)]
 #    col_coef[i, tj:obsi[j]] <- col_coefj
 #   }
 #   tj <- obsi[j] + 1
 #  }
 # }
 # 
 # for (i in 1:nrow(sp_matrix)) {
 #  obsi <- which(sp_matrix[i, ] == 1)
 #  for (j in seq_along(obsi)) {
 #   nyearj <- min(ncol(sp_matrix) - obsi[j] + 1, length(persi_coef))
 #   if(nyearj > 0){
 #    per_coefj <- rev(persi_coef)[1:nyearj]
 #    per_coef[i, obsi[j]:(obsi[j] + nyearj - 1)] <- per_coefj
 #   }
 #  }
 # }
 
 
 
 to_pick <- which(apply(sp_matrix,1,sum)>1)
 detect_data_resp <- c(sp_matrix[to_pick,])
 detect_data_expl <- c(effort_matrix[to_pick,])
 #plot(detect_data_resp ~ detect_data_expl)
 
 m_detect <- glm(detect_data_resp ~ detect_data_expl, family = binomial)
 detect_data_expl <- c(effort_matrix)
 pred_detect <- predict(m_detect,newdata = as.data.frame(detect_data_expl), type = "response")
 pred_detect <- matrix(pred_detect,
                       nrow = nrow(sp_matrix),
                       ncol = ncol(sp_matrix))
 # Assemble global probability matrix
 # NOTE: Several options are provided. Choose the one that suits your needs.
 # Option 1: Adjust using the absence of observations when there is high sampling effort
 # global_proba <- pmax(sp_matrix, col_coef, per_coef, suit_matrix) * (1 - samp_matrix)
 # Option 2: Use the coefficient according to Nicolas
 # global_proba <- pmax(col_coef, per_coef, suit_matrix) * samp_matrix
 # Option 3: Detectability modeled as a function of sampling effort
 global_proba <- pmax(sp_matrix,col_coef, per_coef, suit_matrix) *  (1-pred_detect)
 #plot(global_proba, border=NA)
 
 return(global_proba)
}
#ENLEVER PT_LIST + SP +path au dessus
#sp="sp_1008910"
for( sp in list_species){
 sp.info<- data.frame(species= sp) #put name in df
 random_cara<-sample(1:nrow(strategy), 1)
 sp.info$colonisation_time <- strategy$CT[random_cara]%/%6 #add Colonisation time to df
 sp.info$persistence_time <- strategy$PT[random_cara]%/%6 #add persistence time to df
 saveRDS(sp.info, paste0("saved_data/random_species/", sp, "/table_trend/random_caracteristics.RDS"))

 for (tr in trends){
  trend<- gsub("biais_", "\\1", tr)
  #create a link to the files
  species_files<-list.files(path = paste0("saved_data/random_species/", sp,"/", tr ), pattern = "biais_",full.names = TRUE)
  year_numbers <- as.numeric(gsub(".*year([0-9]+)\\.RDS", "\\1", species_files))
  species_files <- species_files[order(year_numbers)]  # Reorder the files by those extracted numbers

  #read the suitability map
  prob_raster_quantile <- readRDS(paste0("saved_data/random_species/", sp, "/biais_decrease/suitab_random.RDS"))

  #apply the df creation on each time step
  sp_df<-do.call(rbind, lapply(species_files, to_matrix))

  #select just where there is at least one obs
  sp_obs <- sp_df[sp_df$sq_id %in% unique(sp_df$sq_id[sp_df$num_points!=0]), ]
  sp_obs[is.na(sp_obs)] <- 0

  # Convert to matrix format
  sp_matrix_dat <- reshape2::dcast(sp_obs, sq_id ~ year, value.var = "num_points", fill = 0)
  suit_matrix_dat <- reshape2::dcast(sp_obs, sq_id ~ year, value.var = "suitability", fill = 0)
  #samp_matrix_data <- dcast(complete_data, sq_id ~ year, value.var = "proba_coef", fill = 0)
  effort_matrix_dat <- reshape2::dcast(sp_obs, sq_id ~ year, value.var = "sampling_eff", fill = 0)

  # Convert to matrix (excluding sq_id column)
  sp_matrix <- as.matrix(sp_matrix_dat[, -1])
  sp_matrix[sp_matrix > 0] <- 1
  rownames(sp_matrix) <- sp_matrix_dat$sq_id  # Set row names as site IDs

  suit_matrix <- as.matrix(suit_matrix_dat[, -1])
  rownames(suit_matrix) <- suit_matrix_dat$sq_id  # Set row names as site IDs

  # samp_matrix <- as.matrix(samp_matrix_data[, -1])
  # rownames(samp_matrix) <- samp_matrix_data$sq_id  # Set row names as site IDs
  #
  effort_matrix <- as.matrix(effort_matrix_dat[,-1])
  rownames(effort_matrix)

  # plot(sp_matrix,
  # main = paste0("Observations per Site-Year for ",sp, ", ", tr), border = NA)

  #save the three different matrices
  saveRDS(sp_matrix, paste0("saved_data/random_species/", sp, "/table_trend/matrix_sp_", trend,".RDS"))
  saveRDS(suit_matrix, paste0("saved_data/random_species/", sp, "/table_trend/matrix_suit_", trend,".RDS"))
  saveRDS(effort_matrix, paste0("saved_data/random_species/", sp, "/table_trend/matrix_effort_", trend,".RDS"))

  #create the occupancy matrix + save it with the infoflora model
  matrix_sp <- readRDS(paste0("saved_data/random_species/", sp, "/table_trend/matrix_sp_",trend, ".RDS"))
  matrix_effort <- readRDS(paste0("saved_data/random_species/", sp, "/table_trend/matrix_effort_",trend, ".RDS"))
  matrix_suit<- readRDS(paste0("saved_data/random_species/", sp, "/table_trend/matrix_suit_",trend, ".RDS"))
  occupancy_matrix<-occupancy_function2(sp.info,strategy, matrix_sp, matrix_effort,matrix_suit)
  saveRDS(occupancy_matrix, paste0("saved_data/random_species/", sp, "/table_trend/occupancy_matrix_", trend,".RDS"))

  #create the supposed population trend by infoflora model
  occupancy_matrix <- readRDS(paste0("saved_data/random_species/", sp, "/table_trend/occupancy_matrix_", trend, ".RDS"))
  trend_occup<- colSums(occupancy_matrix)
  
  #load the table where the populations and observations per year are stored
  trend_table <- readRDS(paste0("saved_data/random_species/", sp, "/table_trend/", tr, ".RDS"))
  
  #add the trend of infoflora to the table
  trend_table[nrow(trend_table) + 1,] = trend_occup
  row.names(trend_table)[row.names(trend_table) == "5"] <- "Infoflora"
  
  #scale the trends depending of their resolution
  trend_table<-rbind(trend_table,trend_table[1:2,]/trend_table[2,1])
  trend_table<-rbind(trend_table,trend_table[3:4,]/trend_table[4,1])
  trend_table<-rbind(trend_table,trend_table[5,]/trend_table[5,1])
  
  saveRDS(trend_table, paste0("saved_data/random_species/", sp, "/table_trend/table_", trend,".RDS"))
  
  # # Plot the table (Observations)
  # plot(c(1:length(trend_table[1,])), trend_table[3,], type = "o", col = "blue",
  #      xlab = "TimeStep", ylab = "Nb of points",
  #      xlim = c(1, length(trend_table[1,])), ylim = c(0,1.2*max(trend_table[2,])),
  #      main = paste( trend,"species", id,  "over time"), lty= "44")
  # 
  # # Add the second line (Population)
  # lines(c(1:length(trend_table[1,])), trend_table[4,], type = "o", col = "blue")
  # 
  # # Add the sq Observations
  # lines(c(1:length(trend_table[1,])), trend_table[1,], type = "o", col = "orange",lty= "44")
  # 
  # # Add the sq Populations
  # lines(c(1:length(trend_table[1,])), trend_table[2,], type = "o", col = "orange")
  # 
  # lines(c(1:length(trend_table[1,])), trend_table[5,], type = "o", col ="red")
  
  # # Add a legend
  # legend("topright", legend = c("Population", "Observations", "Sq populations", "Sq observations", "Infoflora"),
  #        col = c("blue", "blue", "orange", "orange", "red"), lty= c( "solid", "44", "solid","44", "solid"), pch = 1, ncol =3)
  pdf(paste0("saved_data/random_species/", sp, "/table_trend/graph_", trend,".pdf"), width = 7.4, height = 6.26)
  
  # Plot the table (Observations)
  par(mar=c(5,6,4,1)+.05)
  plot(timescale, trend_table[8,], type = "o", col = "#359B73",
       xlab = "Time", ylab = "Occupied sites [%]",
       xlim = c(1946,2025), ylim = c(0,max(trend_table[6:10,])+0.6),
       main = paste0( trend," species ", sp,  " over time scaled"), lty= "44",  xaxt = "n")
  axis(1, at = timescale, las=2)
  
  # Add the second line (Population)
  lines(timescale, trend_table[9,], type = "o", col = "#359B73")
  
  # Add the sq Observations
  lines(timescale, trend_table[6,], type = "o", col = "#3DB7E9",lty= "44")
  
  # Add the sq Populations
  lines(timescale, trend_table[7,], type = "o", col = "#3DB7E9")
  
  lines(timescale, trend_table[10,], type = "o", col ="#E69F00")
  
  # Add a legend
  legend("topright", legend = c("100m Popu", "100m Obs", "5000m Popu", "5000 Obs", "5000m Infoflora"),
         col = c("#359B73", "#359B73", "#3DB7E9", "#3DB7E9", "#E69F00"), lty= c( "solid", "44", "solid","44", "solid"), pch = 1, ncol =3)
  dev.off()
 }
 # plot(sp_matrix,
 # main = paste0("Observations per Site-Year for ",sp, ", ", trend), border = NA)
 print(paste0(sp, " occupancy + matrices done"))
}
#################################################################################
#multiple models
#################################################################################


list_sp= pt_list
list_sp= list_species
trend= trends
#loop through all the species
mean_caculation<-function(list_sp,trend){
 final_df<-data.frame()
 plot(c(1:14), c(rep(0,14)), type = "o",col= "white",
      xlab = "TimeStep", ylab = "Proportion occupied sites",
      xlim = c(1, 14), ylim = c(0,4),
      main = paste0( "Mean Infoflora models compared to populations ", trend[2], " over Time"), lty= "44")
 for( sp in list_species){
  tr<- gsub("biais_", "\\1", trend[2])
  new_table <- readRDS(paste0("saved_data/random_species/", sp, "/table_trend/table_",tr,".RDS"))
  new_table$species<- sp
  new_table$type<- c("Biais_sq","Population_sq", "Biais", "Population", "Infoflora", "Biais_sq1","Population_sq1", "Biais1", "Population1", "Infoflora1")
  
  # # Add the sq Populations
  lines(c(1:length(new_table[1,])), new_table[7,], type = "l", col = alpha("#E69F00", 0.4), lwd=1.5)
  
  # # Add the second line (Infoflora)
  lines(c(1:length(new_table[1,])), new_table[10,], type = "l", col = alpha("#3DB7E9", 0.4), lwd=1.5)
  
  final_df<- rbind(final_df, new_table[6:10,])
  print(sp)
  
  
 }
 infl_df<-final_df[final_df$type =="Infoflora1",]
 infl_mean <-colSums(infl_df[,1:14])/nrow(infl_df[,1:14])
 
 sq_popu_df<-final_df[final_df$type =="Population_sq1",]
 sq_popu_mean <-colSums(sq_popu_df[,1:14])/nrow(sq_popu_df[,1:14])
 
 lines (infl_mean, type = "l", col = "#2271B2", lwd=2)
 lines(sq_popu_mean, type = "l", col = "#D55E00", lwd=2)
 #Add a legend
 legend("bottomleft", legend = c("Mean 5000m Populations", "5000m virtual species", "mean 5000m Infoflora", "5000m Infoflora model"),
        col = c("#D55E00", "#E69F00", "#2271B2", "#3DB7E9"), bty = "n",lwd=2, cex=0.9, lty=c(1,1,1),ncol =2)
 # legend("bottomleft", legend = c("Mean Sq Populations", " "),
 #         col = c("red", "white"), bty = "n",lwd=2, cex=1.2, lty=c(1,1,1),)
 
}

mean_caculation(list_sp, trends)
