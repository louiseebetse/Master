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
nb_obs<-data.frame(obs_Origin_1994=c(0,0), obs_1995_2000= c(0,0), obs_2001_2006=c(0,0), obs_2007_2012=c(0,0), obs_2013_2018=c(0,0), obs_2019_2024=c(0,0))

for (i in c(1:6)){
 total_sum <- global(my_biais[[i]], fun = "sum", na.rm = TRUE)$sum
 nb_obs[1,i]<-total_sum
}
total_obs<-sum(nb_obs)
for (i in c(1:6)){
 percentage <- nb_obs[1,i]/total_obs
 nb_obs[2,i]<-percentage
}

#################################################################################
#random distributions
#################################################################################
nb_random_species <- 10
for (j in 6:nb_random_species){
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
random.sp4$pa.raster[random.sp4$suitab.raster<0.2]<-0
random.sp4$suitab.raster[random.sp4$suitab.raster<0.05]<-0
#randomly select half of the popu from my random distributions
#new_random_popu<-spatSample(random.sp4$pa.raster, size = global(random.sp4$pa.raster == 1, fun = "sum", na.rm = TRUE)$sum/3 ,
                            #method= "weights", na.rm= TRUE, as.points = TRUE)
#response_realistic <- readRDS("~/Desktop/UNIL/Master/Rstudio/Virtual species/saved_data/random_species/response_realistic.RDS")
#response_realistic$pa.raster<-terra::unwrap(response_realistic$pa.raster)
#response_realistic$suitab.raster<-terra::unwrap(response_realistic$suitab.raster)

random_total_popu<- as.points(random.sp4$pa.raster == 1, values=FALSE)
#saveRDS(random_total_popu, "saved_data/random_species/random_total_popu.RDS")
sample_size<-global(random.sp4$pa.raster == 1, fun = "sum", na.rm = TRUE)$sum%/%3 
#saveRDS(sample_size, "saved_data/random_species/sample_size.RDS" )
suitab_random<-random.sp4$suitab.raster

print(paste("the population nb",j,"has been generated"))

#I decrease the number of sample to see if it's work otherwise it will take too long
#nb_obs<-nb_obs%/%100

#STABLE OBSERVATIONS 90%
biais_stady_sample<- function(repetition, sample_size, sample_vector,proba_raster, biais_raster){
 dir.create(file.path(paste("saved_data/random_species/random_sp_", j,"/biais_stable",sep="")), showWarnings = FALSE)
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
  saveRDS(year,paste("saved_data/random_species/random_sp_", j,"/biais_stable/year",i,".RDS",sep=""))
  
  biais_weights<- terra::extract(biais_raster[[i]], year, xy= TRUE)[, 2]
  #the numbers of observations come from my initial popu(that was divided by three) * the proportions of observations this year
  #then i sample from my population of the year with a weight based on the suitab map
  biais_year<-sample(year, size = as.integer(nb_obs[2,i]*sample_size), prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year,paste("saved_data/random_species/random_sp_", j,"/biais_stable/biais_year",i,".RDS",sep=""))
  
  
 }
 return(invisible(year))
}
biais_stady_sample(6,sample_size,random_total_popu,suitab_random,biais_transformed)

#INCREASE POPULATION WITH OBSERVATIONS 90%
biais_increase_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_increase,biais_raster){
 #set.seed(42)
 dir.create(file.path(paste("saved_data/random_species/random_sp_", j,"/biais_increase",sep="")), showWarnings = FALSE)
 par(mfrow = c(1, 2))
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_weights<-weights
 new_weights[indices_to_keep]<-0
 saveRDS(year,paste("saved_data/random_species/random_sp_", j,"/biais_increase/year1.RDS",sep=""))
 
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = as.integer(nb_obs[2,i]*sample_size), prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year,paste("saved_data/random_species/random_sp_", j,"/biais_increase/biais_year1.RDS",sep=""))
 
 
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
  saveRDS(year,paste("saved_data/random_species/random_sp_", j,"/biais_increase/year",i+1,".RDS",sep=""))
  
  
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size+i*num_to_add, prob = biais_weights, replace = TRUE)
  biais_year<-sample(year, size = as.integer(nb_obs[2,i+1]*length(year)), prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year,paste("saved_data/random_species/random_sp_", j,"/biais_increase/biais_year",i+1,".RDS",sep=""))
  
 }
 return(invisible(year))
}
biais_increase_sample(6,sample_size,random_total_popu,suitab_random,50,biais_transformed)

#DECREASE POPULATION WITH OBSERVATIONS 90%
biais_decrease_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_reduce,biais_raster){
 #set.seed(42)
 dir.create(file.path(paste("saved_data/random_species/random_sp_", j,"/biais_decrease",sep="")), showWarnings = FALSE)
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 inverse_weights <- 1 - weights
 
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_inverse_weights<-inverse_weights[indices_to_keep]
 saveRDS(year,paste("saved_data/random_species/random_sp_", j,"/biais_decrease/year1.RDS",sep=""))
 
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = as.integer(nb_obs[2,i]*sample_size), prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year,paste("saved_data/random_species/random_sp_", j,"/biais_decrease/biais_year1.RDS",sep=""))
 
 
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
  saveRDS(year,paste("saved_data/random_species/random_sp_", j,"/biais_decrease/year",i+1,".RDS",sep=""))
  
  
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size-i*num_to_delete, prob = biais_weights, replace = TRUE)
  biais_year<-sample(year, size = as.integer(nb_obs[2,i+1]*length(year)), prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year,paste("saved_data/random_species/random_sp_", j,"/biais_decrease/biais_year",i+1,".RDS",sep=""))

 }
 return(invisible(year))
}
biais_decrease_sample(6,sample_size,random_total_popu,suitab_random,50,biais_transformed)

print(paste("stable + increase + decrease of random sp",j,"has been generated"))

saveRDS(suitab_random, paste("saved_data/random_species/random_sp_", j,"/biais_stable/suitab_random.RDS",sep=""))
saveRDS(suitab_random, paste("saved_data/random_species/random_sp_", j,"/biais_increase/suitab_random.RDS",sep=""))
saveRDS(suitab_random, paste("saved_data/random_species/random_sp_", j,"/biais_decrease/suitab_random.RDS",sep=""))
}

#### table of 5km square for popu and observations
#acess the data saved
# decrease_table<- function(raster5x5 = my_biais$obs_Origin_1994, file.path = "saved_data/random_species/random_sp_2/biais_increase" ){
#  all_the_files<-list.files(path = file.path,
#                            pattern = NULL,
#                            full.names = TRUE)
#  data <- matrix(NA, nrow = 4, ncol = length(all_the_files)/2) #create the matrix
#  colnames(data) <- paste0("TimeStep", 1:6)
#  rownames(data) <- c("Biais_sq", "Population_sq", "Biais", "Population")
#  data <- as.data.frame(data)
#  for (i in 1:(length(all_the_files)/2)){
#   biais_year<-readRDS(all_the_files[i])
#   year<-readRDS(all_the_files[i+length(all_the_files)/2+1])
# 
#   biais_cell_ids <- cellFromXY(raster5x5, crds(biais_year))
#   cell_ids <- cellFromXY(raster5x5, crds(year))
# 
#   # Count the number of observations per cell
#   counts_cell <- length(table(cell_ids))
#   counts_biais <- length(table(biais_cell_ids))
# 
#   data[1, i] <- counts_biais #nb of square of observations
#   data[2, i] <- counts_cell #nb of square of populations
#   data[3, i] <- length(biais_year)# nb of observations
#   data[4, i] <- length(year)# nb of populations
#  }
#  return(data)
# }
# decrease_table(file.path = "saved_data/random_species/random_sp_2/biais_increase")
