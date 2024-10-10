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
#rasterise my variables
my_variables<- rast(my_variables)
cat(crs(my_variables))
#change the names of my variables
names(my_variables) = c("ai","bio4", "bio6", "bio15","gdd3","ph", "pop2000", "pop5000", "rad", "dem", "slope", "canopy")
#compute the curvature
my_variables$dem<-terrain(my_variables$dem, v= "TPI")
##plot a raster
#plot(my_variables)

#export the swiss shapefile
swiss_shape<-st_read("swiss_shapefile.gpkg")
swiss_shape <- st_transform(swiss_shape, crs = 2056)

#export my species data from my file
load("SPImaster/sp/1008910")

## CLEAR THE DATA
C.kitaibelii<-clear(my.sp) #clear my species data using clear()
C.kitaibelii <- SpatialPointsDataFrame(coords=C.kitaibelii[,c("x", "y")], data=C.kitaibelii[ , -c(6, 7)]) #creating a SPDF with just the coordinates
C.kitaibelii<- remove.duplicates(C.kitaibelii, zero = 300) #remove duplicate presence in a 300 radius?
C.kitaibelii<- as.data.frame(C.kitaibelii) #come back to a data frame
coord<- C.kitaibelii[,c("x", "y")]
#plot(C.kitaibelii)


#extract the variables coord at the species presence
variablesValues<-data.frame(extract(my_variables,coord))
head(variablesValues)

#put the species data in the same table
spData<- na.omit(data.frame(coord, variablesValues))
head(spData)
spData<-spData[, !names(spData) %in% "ID"]#remove the id column
spData$presence <- rep(1, times = length(spData$x))#put a presence line

#generate background points 
bg_rand <- terra::spatSample(my_variables$ai, 10000, "random", na.rm=T, as.points=TRUE)
#plot(swiss_shape$geom,col='lightgrey')
#points(coord,pch='+',col='red')# one way to plot presences
#points(C.kitaibelii$x, C.kitaibelii$y, col='red',  pch=19) #two ways to plot presences
#points(bg_rand,pch=13,cex=0.3)

#create a dataframe with presences and pseudo absences
sp_env <- data.frame(coord, occ=1)# create dataframe with presence
bg_rand_buf_df <- data.frame(terra::geom(bg_rand)[,c('x','y')])# create dataframe with absence
bg_rand_buf_df$occ <- 0 #add the 0 to the dataframe with absences
sp_env <- rbind(sp_env, bg_rand_buf_df) #bind the two dataframe in one big 
sp_env <- cbind(sp_env, terra::extract(x = my_variables , y = sp_env[,c('x','y')], cells=T) )
sp_env<-sp_env[, !names(sp_env) %in% c("ID","cell", "pop2000", "pop5000", "dem")]
sp_env <- na.omit(sp_env)# remove the NA of the data


#test<- sp_env[c(1:200,900:1100),]
#args<- cov.sel(T=sp_env$occ, Y= sp_env$x , X= sp_env[,4:14], type = "np", alpha = 0.3)
#glmModAIC<-stepAIC(glmStart, glm.formula, direction= "both", trace= FALSE, k=2)
#summary(args)

#create the glm
glm.formula<- bm_MakeFormula("sp_env$occ",sp_env[, c("bio6", "canopy", "rad", "ph", "ai", "bio4", "bio15", "gdd3","slope")], "quadratic", interaction.level= 0)
glm1<- glm(glm.formula, data=sp_env, family= "binomial")


#create formula for virtual species manually
var.names<-c("bio6", "canopy", "rad", "ph", "ai", "bio4", "bio15", "gdd3","slope")
inter= glm1$coefficients[1]#create the interception to y
medians<- apply(sp_env[, var.names], 2, median)#create the medians of all variables
medians_table<- data.frame(sapply(medians,function(x)rep(x,100)))
median<- medians_table[1,]#have a data frame with one line
coef<-glm1$coefficients# create a dataframe with just the coef (a in ax +b) from the glm
#function to have h (intercept + coef*median in every other variable except the first one)

#function to create the different 
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
  funcs[[i]] <- function(x,h,b,c) {
   h + b*x + c*x^2
  }
 }
 df<- data.frame(variable= var_names, h= H, b= B, c=C, stringsAsFactors = FALSE)
 df$functions <- funcs
 row.names(df)<-var_names
 return(df)
}

fun1<- df$functions[[1]]
fun2<- df$functions[[2]]
fun3<- df$functions[[3]]
fun4<- df$functions[[4]]
fun5<- df$functions[[5]]
fun6<- df$functions[[6]]
fun7<- df$functions[[7]]
fun8<- df$functions[[8]]
fun9<- df$functions[[9]]

response.variable<- function(median, coef, interest_variable){
 h= as.numeric(coef[1])#h=intercept
 median<- median[, !names(median) %in% interest_variable]# remove bio6 from the medians
 coef<-coef[!names(coef) %in% c("(Intercept)", interest_variable, paste0("I(", interest_variable, "^2)"))]# remove intercept and bio6 coef
 for (i in 1:length(pro_median)){
  h=as.numeric(h+ coef[i*2-1]*median[i]+ coef[i*2]*median[i]^2)
 }
 return (h)
}

h= response.variable(median, coef, "bio6")
b=as.numeric(coef[2])
c=as.numeric(coef[3])
var1.function <- function(x, h,b,c){
 h + b*x + c*x^2
}

response.variable2<- function(median, coef){
 h= as.numeric(coef[1])
 median<- median[,-2]
 interest_coef<-coef[4:5]
 coef<-coef[-c(1,4,5)]
 for (i in 1:length(median)){
  h=as.numeric(h+ coef[i*2-1]*median[i]+ coef[i*2]*median[i]^2)
 }
 return (h)
}
h2= response.variable2(median, coef)
b2=as.numeric(coef[4])
c2=as.numeric(coef[5])
var2.function <- function(x, h,b,c){
 h + b*x + c*x^2
}


my.responses <- formatFunctions(bio6 = c(fun = 'test1', h=df$h[1],b=df$b[1], c=df$c[1]),
                                canopy = c(fun = 'test2', h=df$h[2],b=df$b[2], c=df$c[2]),
                                rad = c(fun = 'test2', h=df$h[3],b=df$b[3], c=df$c[3]),
                                ph = c(fun = 'test2', h=df$h[4],b=df$b[4], c=df$c[4]),
                                ai = c(fun = 'test2', h=df$h[5],b=df$b[5], c=df$c[5]),
                                bio4 = c(fun = 'test2', h=df$h[6],b=df$b[6], c=df$c[6]),
                                bio15 = c(fun = 'test2', h=df$h[7],b=df$b[7], c=df$c[]),
                                gdd3 = c(fun = 'test2', h=df$h[8],b=df$b[8], c=df$c[1]),
                                slope = c(fun = 'test2', h=df$h[9],b=df$b[9], c=df$c[1]),)

generateSpFromFun(raster.stack = my_variables[[c("bio6", "canopy")]],
                  parameters = my.responses, plot = TRUE)


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
