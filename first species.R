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

#generate background points 
bg_rand <- terra::spatSample(my_variables$ai, 10000, "random", na.rm=T, as.points=TRUE)
plot(swiss_shape$geom,col='lightgrey')
points(coord,pch='+',col='red')# one way to plot presences
points(C.kitaibelii$x, C.kitaibelii$y, col='red',  pch=19) #two ways to plot presences
points(bg_rand,pch=13,cex=0.3)

#create a dataframe with presences and pseudo absences
sp_env <- data.frame(coord, occ=1)# create dataframe with presence
bg_rand_buf_df <- data.frame(terra::geom(bg_rand)[,c('x','y')])# create dataframe with absence
bg_rand_buf_df$occ <- 0 #add the 0 to the dataframe with absences
sp_env <- rbind(sp_env, bg_rand_buf_df) #bind the two dataframe in one big 
sp_env <- cbind(sp_env, terra::extract(x = my_variables , y = sp_env[,c('x','y')], cells=T) )

#sp_env<-sp_env[,-9]
sp_env <- na.omit(sp_env)# remove the NA of the data
#test<- sp_env[c(1:200,900:1100),]
#args<- cov.sel(T=sp_env$occ, Y= sp_env$x , X= sp_env[,4:14], type = "np", alpha = 0.3)
#glmModAIC<-stepAIC(glmStart, glm.formula, direction= "both", trace= FALSE, k=2)
#summary(args)

#create the glm
glm.formula<- bm_MakeFormula("sp_env$occ",sp_env[, c("bio6", "canopy", "rad", "dem", "ai", "bio4", "bio15", "gdd3","slope", "cell")], "quadratic", interaction.level= 0)
glm1<- glm(glm.formula, data=sp_env, family= "binomial")


# out1.glm= NULL
# var.names<- c("bio6", "canopy", "rad", "dem", "ai", "bio4", "bio15", "gdd3","slope", "cell")
# medians<- apply(sp_env[, var.names], 2, median)
# medians_table<- data.frame(sapply(medians,function(x)rep(x,100)))
# 
# 
# for(i in 1:10){
#  foc.var <- sp_env[,var.names[i]]
#  new.data <- medians_table
#  var.new <- seq(min(foc.var), max(foc.var), length=100)
#  new.data[,i] <- var.new
#  pred.glm <- predict(glm1, newdata=new.data, type="response")
#  #glm2 <- predict( glmStart, newdata=new.data, type="response")
#  #pred.glob<- c(pred.glm,glm2)
#  tmp1 <- cbind(Occ.prob=pred.glm,Env.val=var.new)
#  tmp2 <- data.frame(cbind(Algorithm=rep("GLM1",each=100),
#                           Var.name=var.names[i]))
#  out1.glm <- rbind(out1.glm,cbind(tmp1,tmp2))
# }
# #plot the predictions
# resp1.glm <- ggplot(out1.glm, aes(x=Env.val, y=Occ.prob, color=Algorithm,
#                                   linetype=Algorithm)) +
#  geom_line(lwd=1) +
#  scale_color_manual(values=c("#007991", "#439A86", "#BCD8C1"))+
#  facet_wrap(~Var.name,ncol=2,scale="free_x")+
#  theme_bw()+
#  labs(x="Value",y="Occurence probability")+
#  theme_bw()
# print(resp1.glm)
# 
# Pred_test <- predict(glmModAIC, calib, type="response")
# plot(Pred_test)

#create formula for virtual species manually
inter= glm1$coefficients[1]#create the interception to y
medians<- apply(sp_env[, var.names], 2, median)#create the medians of all variables
medians_table<- data.frame(sapply(medians,function(x)rep(x,100)))
median<- medians_table[1,]#have a data frame with one line
coef<-glm1$coefficients# create a dataframe with just the coef (a in ax +b) from the glm
#function to have h (intercept + coef*median in every other variable except the first one)
response.variable<- function(median, coef){
 h= as.numeric(coef[1])#h=intercept
 median<- median[,-1]# remove bio6 from the medians
 coef<-coef[4:length(coef)]# remove intercept and bio6 coef
 for (i in 1:length(median)){
  h=as.numeric(h+ coef[i*2-1]*median[i]+ coef[i*2]*median[i]^2)
  #print(h)
 }
 #response_variable<- as.numeric(h + interest_coef[1]*x +interest_coef[2]*x)
 return (h)
}
h= response.variable(median, coef)
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


my.responses <- formatFunctions(bio6 = c(fun = "var1.function", h=h,b=b, c=c),
                                canopy = c(fun = "var2.function", h=h2,b=b2, c=c2))

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
