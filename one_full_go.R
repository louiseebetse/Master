library(virtualspecies)
library(terra)
library(geodata)
library(git2r)
library(sf)
library(gtools) # V.3.9.4
library(data.table)
library(sp)
library(biomod2)
library(MASS)

#set the session to the right place
setwd("/Users/louise/Desktop/UNIL/Master/Rstudio/Virtual species")
dir.create(file.path("saved_data/suitab_map"), showWarnings = FALSE)

#export my variables from my file
my_variables <- mixedsort(list.files(path = "variables",
                                     pattern = NULL,
                                     full.names = TRUE))
my_variables<- rast(my_variables)
#cat(crs(my_variables))
names(my_variables) = c("ai","bio4", "bio6", "bio15","gdd3", "ph", "pop2000", "pop5000", "rad", "dem", "slope", "canopy")
my_variables$rad <- project(my_variables$rad, "EPSG:2056")
my_variables$dem<-terrain(my_variables$dem, v= "TPI")
my_variables<- subset(my_variables, c("ai","bio4", "bio6", "bio15","gdd3", "ph", "dem", "slope", "canopy")) #dem, slope canopy
my_variables<- aggregate(my_variables, fact= 4)#aggregate to have a 100m resolution

#export the swiss shapefile
swiss_shape<-st_read("swiss_shapefile.gpkg")
swiss_shape <- st_transform(swiss_shape, crs = 2056)

#export my species data from my file
#load("SPImaster/sp/1008880")
#sp<-my.sp
#load("SPImaster/sp/1008910")
# Initialize a list to store the data
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

pt_list<- data_list[c(2,17)]
#name<- names(pt_list)
name<- names(data_list)

fun1<- function(x,h,b,c) {
 h + b*x + c*x^2
}

#load the file with the biais raster
my_biais <- mixedsort(list.files(path = "biais",
                                 pattern = NULL,
                                 full.names = TRUE))
my_biais<- rast(my_biais)
my_biais<-terra::unwrap(my_biais)
#head(my_biais)
names(my_biais)<-c("obs_Origin_1994", "obs_1995_2000", "obs_2001_2006", "obs_2007_2012",
                   "obs_2013_2018", "obs_2019_2024",
                   "Nobsvs_Origin_1994", "Nobsvs_1995_2000", "Nobsvs_2001_2006", "Nobsvs_2007_2012",
                   "Nobsvs_2013_2018", "Nobsvs_2019_2024",
                   "prop_pot_Origin-1994", "prop_pot_1995-2000", "prop_pot_2001-2006",  
                   "prop_pot_2007-2012", "prop_pot_2013-2018", "prop_pot_2019-2024")


#apply the log scale to my biais raster
biais_transformed<- my_biais
for (i in 1:length(names(my_biais))){
 biais_transformed[[i]]<-log(my_biais[[i]] + 1)
}

#biais_stady_sample(repetition=6,size_popus[sp],sample2,prob_raster_quantile,biais_transformed,path_1)


#functions to increase, decrease + stady:
biais_stady_sample<- function(repetition, sample_size, sample_vector,proba_raster, biais_raster,path){
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
  biais_year<-sample(year, size = date_vector[i], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste(path, "/biais_year",i,".RDS", sep=""))#save each year in a file
  
  
 }
 return(invisible(year))
}
biais_increase_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_increase,biais_raster,path,obs_by_year){
 #set.seed(42)
 par(mfrow = c(1, 2))
 weights <- terra::extract(proba_raster, sample_vector, xy= TRUE)[, 2]# i extract the weights for my sample from my proba raster
 indices_to_keep<-sample(nrow(sample_vector), size = sample_size, 
                         prob = weights, replace = FALSE)
 year<-sample_vector[indices_to_keep,]#recreate a spatvector from the desired individuals
 new_weights<-weights
 new_weights[indices_to_keep]<-0
 saveRDS(year, paste(path, "/year1.RDS", sep="")) #save each year in a file
 
 biais_weights<- terra::extract(biais_raster[[1]], year, xy= TRUE)[, 2]
 biais_year<-sample(year, size = date_vector[1], prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year, paste(path, "/biais_year1.RDS", sep=""))
 
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
  saveRDS(year, paste(path, "/year",i+1,".RDS", sep=""))#save each year in a file
  
  biais_weights<- terra::extract(biais_raster[[i+1]], year, xy= TRUE)[, 2]
  #biais_year<-sample(year, size = sample_size+i*num_to_add, prob = biais_weights, replace = TRUE)
  biais_year<-sample(year, size = date_vector[i+1], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste(path, "/biais_year",i+1,".RDS", sep=""))#save each year in a file
 }
 return(invisible(year))
}
biais_decrease_sample<- function(repetition, sample_size,sample_vector,proba_raster, max_reduce,biais_raster,path,obs_by_year){
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
 biais_year<-sample(year, size = date_vector[1], prob = biais_weights, replace = TRUE)
 unique_biais_year <- terra::unique(biais_year)
 saveRDS(unique_biais_year, paste(path, "/biais_year1.RDS", sep=""))
 num_to_delete <- (max_reduce/(repetition-1))/100*sample_size
 if (num_to_delete <1) {num_to_delete<-1}
 
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
  biais_year<-sample(year, size = date_vector[1+i], prob = biais_weights, replace = TRUE)
  unique_biais_year <- terra::unique(biais_year)
  par(mfg = c(1, 2))
  points(unique_biais_year, col= col[i+1], pch= 20, cex =0.5)#make points on the map
  saveRDS(unique_biais_year, paste(path, "/biais_year",i+1,".RDS", sep=""))#save each year in a file
  
 }
 return(invisible(year)) 
}

#remove useless variables
rm(data,my.sp,file_name,file,dwld.date)

my_distribution<- function(list_presences){
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
  save(coord,file = paste("saved_data/species_1008910.RData"))
  
  #create a dataframe with presences and pseudo absences
  bg_rand <- terra::spatSample(my_variables$dem, 10000, "random", na.rm=T, as.points=TRUE)
  sp_env <- data.frame(coord, occ=1)# create dataframe with presence
  bg_rand_buf_df <- data.frame(terra::geom(bg_rand)[,c('x','y')])# create dataframe with absence
  bg_rand_buf_df$occ <- 0 #add the 0 to the dataframe with absences
  sp_env <- rbind(sp_env, bg_rand_buf_df) #bind the two dataframe in one big 
  sp_env <- cbind(sp_env, terra::extract(x = my_variables , y = sp_env[,c('x','y')], cells=T) )
  sp_env<-sp_env[, !names(sp_env) %in% c("ID","cell", "pop2000", "pop5000", "rad")]
  sp_env <- na.omit(sp_env)# remove the NA of the data
  
  #create the glm
  glm_formula<- sp_env$occ ~ 1 + bio6 + I(bio6^2) + canopy + I(canopy^2) + dem + 
   I(dem^2) + ph + I(ph^2) + ai + I(ai^2) + bio4 + I(bio4^2) + 
   bio15 + I(bio15^2) + gdd3 + I(gdd3^2) + slope + I(slope^2)
  glm1<- glm(glm_formula, data=sp_env, family= "binomial")
  #sp_prediction<-predict.glm(glm1,newdata= my_variables, type= "response")
  
  #create a vector with the names of the variables
  var.names<-c("bio6", "canopy", "dem", "ph", "ai", "bio4", "bio15", "gdd3","slope")
  
  #create thedataframe with the medians
  medians <- data.frame(t(apply(sp_env[, var.names], 2, median)))
  
  
  # create a dataframe with just the coef (a in ax +b) from the glm
  coef<-glm1$coefficients
  
  #a function that create a dataframe to have the h,b,c and the function for each variable
  function_creator <- function(median, coef, var_names){
   funcs <- vector("list", length(var_names))
   H<-c()
   B<-c()
   C<-c()
   for (i in 1:length(var_names)) {
    interest_variable = var_names[i]
    h= as.numeric(coef[1])
    pro_median<- median[, !names(median) %in% interest_variable]# remove bio6 from the medians
    pro_coef<-coef[!names(coef) %in% c("(Intercept)", interest_variable, paste0("I(", interest_variable, "^2)"))]# remove intercept and bio6 coef
    for (j in 1:length(pro_median)){
     h=as.numeric(h+ pro_coef[j*2-1]*median[j]+ pro_coef[j*2]*median[j]^2)
    }
    H<- c(H, h)
    B=c(B, as.numeric(coef[2*i]))
    C=c(C, as.numeric(coef[2*i+1]))
   }
   df<- data.frame(variable= var_names, h= H, b= B, c=C, stringsAsFactors = FALSE)
   row.names(df)<-var_names
   return(df)
  }
  df<- function_creator(medians, coef, var.names)
  
  my.responses <- formatFunctions(bio6 = c(fun = 'fun1', h=df$h[1],b=df$b[1], c=df$c[1]),
                                  canopy = c(fun = 'fun1', h=df$h[2],b=df$b[2], c=df$c[2]),
                                  dem = c(fun = 'fun1', h=df$h[3],b=df$b[3], c=df$c[3]),
                                  ph = c(fun = 'fun1', h=df$h[4],b=df$b[4], c=df$c[4]),
                                  ai = c(fun = 'fun1', h=df$h[5],b=df$b[5], c=df$c[5]),
                                  bio4 = c(fun = 'fun1', h=df$h[6],b=df$b[6], c=df$c[6]),
                                  bio15 = c(fun = 'fun1', h=df$h[7],b=df$b[7], c=df$c[7]),
                                  gdd3 = c(fun = 'fun1', h=df$h[8],b=df$b[8], c=df$c[8]),
                                  slope = c(fun = 'fun1', h=df$h[9],b=df$b[9], c=df$c[9]),)
  map_probability<- generateSpFromFun(raster.stack = my_variables[[c("bio6", "canopy","dem","ph", "ai", "bio4", "bio15", "gdd3", "slope")]],
                                      parameters = my.responses, plot = TRUE)
  
  saveRDS(map_probability,file = paste("saved_data/suitab_map/suitab_map_",name[sp],".RDS", sep=""))
  print(paste("the suitab map of",name[sp],"is done"))
 #}
 #return (size_popus)
#}


#size_popus<-my_distribution(pt_list)

# for(sp in list.files(path = "saved_data/suitab_map",
#                      pattern = NULL,
#                      full.names = TRUE)){
############################################################################
#Sample real popu
############################################################################
 # i=1
 # RDSproba <- readRDS(sp)
 # RDSproba[["suitab.raster"]]<-terra::unwrap(RDSproba[["suitab.raster"]])
 # terra::plot(RDSproba[["suitab.raster"]])
 # 
 ##extract probability for presences of the real species
 #prob_real_species <- terra::extract(RDSproba[["suitab.raster"]], coord, ID= F)
 prob_real_species <- terra::extract(map_probability[["suitab.raster"]], coord, ID= F)
 ##find the probability value under above which we find x% of the sample
 quantile<- quantile(prob_real_species$`VSP suitability`, 0.3, na.rm= T)
 ##i put all the raster values below this quantile values at 0
 prob_raster_quantile<-map_probability[["suitab.raster"]]
 prob_raster_quantile[prob_raster_quantile < quantile] <- 0

 
 ##sample presences from the modified raster(prob_raster_quantile) using 
 ##weights of probability from the raster output spatvector
 sample2<- spatSample(prob_raster_quantile,size= size_popus[sp]*3, method= "weights", na.rm= TRUE, as.points = TRUE)
 print(paste("the prob_raster_quantile and the sample of",name[sp],"is done"))
 
############################################################################
#Sample with biais
############################################################################
#load the species 
#load(paste("SPImaster/sp/",name[i], sep=""))

# Convert the "date" column to Date format if it's not already
 list_presences[[sp]]$date <- as.Date(list_presences[[sp]]$date)
# Define the date breaks
 breaks <- as.Date(c("1994-12-31", "2000-12-31", "2006-12-31", 
                    "2012-12-31", "2018-12-31", "2024-12-31"))
# Use cut() to classify the dates into groups
 categories <- cut(list_presences[[sp]]$date, breaks = c(as.Date("1900-01-01"), breaks), 
                  labels = c("Before 1994", "1995-2000", "2001-2006", 
                             "2007-2012", "2013-2018", "2019-2024"), right = TRUE)

# Count the number of observations in each category & Convert to a numeric vector
 date_vector <- as.numeric(table(categories))
#I now have a vector with the number of observations of my plant per time period
#print(date_vector)
 print(paste("the date vector of",name[sp],"is done"))

 

 #plot(my_biais[[7]])
 dir.create(file.path(paste("saved_data/species/sp_",name[sp],sep="")), showWarnings = FALSE)
 
 #Stady popu
 path_1= paste("saved_data/species/sp_",name[sp],"/biais_stable", sep="")
 dir.create(file.path(path_1), showWarnings = FALSE)
 biais_stady_sample(repetition=6,size_popus[sp],sample2,prob_raster_quantile,biais_transformed,path_1)

 #increase popu
 path_1= paste("saved_data/species/sp_",name[sp],"/biais_increase", sep="")
 dir.create(file.path(path_1), showWarnings = FALSE)
 biais_increase_sample(6,size_popus[sp],sample2,prob_raster_quantile,50,biais_transformed,path_1, date_vector)
 
 #decrease popu
 path_1= paste("saved_data/species/sp_",name[sp],"/biais_decrease", sep="")
 dir.create(file.path(path_1), showWarnings = FALSE)
 biais_decrease_sample(6,size_popus[sp],sample2,prob_raster_quantile, 50, biais_transformed, path_1, date_vector)
 
 #low increase popu
 path_1= paste("saved_data/species/sp_",name[sp],"/low_biais_increase", sep="")
 dir.create(file.path(path_1), showWarnings = FALSE)
 biais_increase_sample(6,size_popus[sp],sample2,prob_raster_quantile, 10, biais_transformed,path_1,date_vector)
 
 #low decrease popu
 path_1= paste("saved_data/species/sp_",name[sp],"/low_biais_decrease", sep="")
 dir.create(file.path(path_1), showWarnings = FALSE)
 biais_decrease_sample(6,size_popus[sp],sample2,prob_raster_quantile, 10, biais_transformed,path_1, date_vector)
 
 print(paste("stable + increase + decrease of", name[sp], "is done"))

 saveRDS(prob_raster_quantile, paste("saved_data/species/sp_",name[sp],"/biais_stable/prob_raster_quantile.RDS",sep=""))
 saveRDS(prob_raster_quantile, paste("saved_data/species/sp_",name[sp],"/biais_increase/prob_raster_quantile.RDS",sep=""))
 saveRDS(prob_raster_quantile, paste("saved_data/species/sp_",name[sp],"/biais_decrease/prob_raster_quantile.RDS",sep=""))
 saveRDS(prob_raster_quantile, paste("saved_data/species/sp_",name[sp],"/low_biais_increase/prob_raster_quantile.RDS",sep=""))
 saveRDS(prob_raster_quantile, paste("saved_data/species/sp_",name[sp],"/low_biais_decrease/prob_raster_quantile.RDS",sep=""))

}
}
my_distribution(data_list)
#my_distribution(pt_list)

#squared distribution

as.square<- function(popu,obs){
 #### put back in a 5km raster 
 # Extract cell IDs for each observation
 cell_ids <- cellFromXY(my_biais$obs_Origin_1994, crds(popu)) 
 biais_cell_ids <- cellFromXY(my_biais$obs_Origin_1994, crds(obs))
 #obs_Origin n'a pas besoin d'etre changer d'une iteration à l'autre
 
 # Count the number of observations per cell
 counts_cell <- table(cell_ids)
 counts_biais <- table(biais_cell_ids)
 #ICI on sait deja la taille
 
 # Create an empty raster with the same structure as your original raster
 counts_raster <- my_biais$obs_Origin_1994
 counts_raster_biais <- my_biais$obs_Origin_1994
 values(counts_raster) <-ifelse(is.na(values(counts_raster)), NA, 0)  # Initialize all cells with 0
 values(counts_raster_biais) <- ifelse(is.na(values(counts_raster)), NA, 0) # Initialize all cells with 0
 
 # Assign counts to the appropriate cells
 values(counts_raster)[as.numeric(names(counts_cell))] <- counts_cell
 values(counts_raster_biais)[as.numeric(names(counts_biais))] <- counts_biais
 
 # Plot the raster of counts
 #plot(counts_raster, main = "Population Counts per 5km Pixel")
 #plot(counts_raster_biais, main = "Observation Counts per 5km Pixel")
 
 #test[[1]]= counts_raster
 #test[[2]]= counts_raster_biais
 return( c(counts_raster, counts_raster_biais))
 
}

species_files<-list.files(path = "saved_data/species", pattern = NULL,full.names = TRUE)
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
  saveRDS(squares[[1]], paste(species_files[sp],"/sq_biais_increase/year",yr,".RDS", sep=""))
  saveRDS(squares[[2]], paste(species_files[sp],"/sq_biais_increase/biais_year",yr,".RDS", sep=""))
  
  #decrease
  popu_decrease<- readRDS(paste(species_files[sp],"/biais_decrease/year", yr, ".RDS", sep=""))
  obs_decrease<-readRDS(paste(species_files[sp],"/biais_decrease/biais_year", yr, ".RDS", sep=""))
  squares<- as.square(popu_decrease,obs_decrease)
  saveRDS(squares[[1]], paste(species_files[sp],"/sq_biais_decrease/year",yr,".RDS", sep=""))
  saveRDS(squares[[2]], paste(species_files[sp],"/sq_biais_decrease/biais_year",yr,".RDS", sep=""))
  
 }
} 
