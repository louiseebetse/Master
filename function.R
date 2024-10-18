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
library(biomod2)
library(MASS)

#set the session to the right place
setwd("/Users/louise/Desktop/UNIL/Master/Rstudio/Virtual species")

#export my variables from my file
my_variables <- mixedsort(list.files(path = "variables",
                                     pattern = NULL,
                                     full.names = TRUE))
my_variables<- rast(my_variables)
cat(crs(my_variables))
names(my_variables) = c("ai","bio4", "bio6", "bio15","gdd3", "ph", "pop2000", "pop5000", "rad", "dem", "slope", "canopy")
my_variables$dem<-terrain(my_variables$dem, v= "TPI")
my_variables<- subset(my_variables, c("ai","bio4", "bio6", "bio15","gdd3", "ph", "dem", "slope", "canopy")) 


#export the swiss shapefile
swiss_shape<-st_read("swiss_shapefile.gpkg")
swiss_shape <- st_transform(swiss_shape, crs = 2056)

#export my species data from my file
#load("SPImaster/sp/1008880")

folder_path<- "SPImaster/sp"

file_list <- list.files(path = folder_path, full.names = TRUE)

# Initialize a list to store the data
data_list <- list()

# Loop through each file
for (file in file_list) {
 # Read the file
 load(file)
 data<-my.sp
 file_name <- tools::file_path_sans_ext(basename(file))
 # Store the data in the list with the file name as the key
 data_list[[file_name]] <- data
}


#pt_list<- data_list[c(2,17)]


fun1<- function(x,h,b,c) {
 h + b*x + c*x^2
}

my_distribution<- function(list_presences){
 list_distribution<- list()
 name<- names(pt_list)
 for (sp in seq_along(pt_list)) {
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
  
  species<-clear(pt_list[[sp]]) #clear my species data using clear()
  species <- SpatialPointsDataFrame(coords=species[,c("x", "y")], data=species[ , -c(6, 7)]) #creating a SPDF with just the coordinates
  species<- remove.duplicates(species, zero = 1000) #remove duplicate presence in a 300 radius?
  species<- as.data.frame(species) #come back to a data frame
  coord<- species[,c("x", "y")]
  
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
  
  #create the functions for the later use
  #fun1<- function(x,h,b,c) {
  # h + b*x + c*x^2
  #}
  my.responses <- formatFunctions(bio6 = c(fun = 'fun1', h=df$h[1],b=df$b[1], c=df$c[1]),
                                  canopy = c(fun = 'fun1', h=df$h[2],b=df$b[2], c=df$c[2]),
                                  dem = c(fun = 'fun1', h=df$h[3],b=df$b[3], c=df$c[3]),
                                  ph = c(fun = 'fun1', h=df$h[4],b=df$b[4], c=df$c[4]),
                                  ai = c(fun = 'fun1', h=df$h[5],b=df$b[5], c=df$c[5]),
                                  bio4 = c(fun = 'fun1', h=df$h[6],b=df$b[6], c=df$c[6]),
                                  bio15 = c(fun = 'fun1', h=df$h[7],b=df$b[7], c=df$c[7]),
                                  gdd3 = c(fun = 'fun1', h=df$h[8],b=df$b[8], c=df$c[8]),
                                  slope = c(fun = 'fun1', h=df$h[9],b=df$b[9], c=df$c[9]),)
  
  map_distribution<- generateSpFromFun(raster.stack = my_variables[[c("bio6", "canopy","dem","ph", "ai", "bio4", "bio15", "gdd3", "slope")]],
                    parameters = my.responses, plot = TRUE)
  list_distribution[[name[sp]]]<-map_distribution
  print(name[sp])
 }
 return(list_distribution)
 }

my_dist_species<-my_distribution(data_list)