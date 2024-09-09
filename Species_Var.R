library(virtualspecies)
library(terra)
library(geodata)
library(git2r)
library(sf)
library(gtools) # V.3.9.4

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
plot(my_variables$dem)
#my_variables

#create the example parameters
my_parameters <- formatFunctions(ai = c(fun= "dnorm", mean = 1500, sd= 2000),
                                 bio4 =c(fun= "dnorm", mean= 57500,sd = 7500),
                                 bio6 =c(fun= "logisticFun", beta = 50, alpha= -400),
                                 bio15 = c(fun= "dnorm", mean = 2000, sd= 2000),
                                 gdd3 =c(fun= "logisticFun", beta = 500, alpha= 7000),
                                 pop2000= c(fun= "logisticFun", beta = -100, alpha = 900000),
                                 pop5000= c(fun= "logisticFun", beta = -100, alpha = 2500000),
                                 rad= c(fun= "dnorm", mean = 50, sd= 20),
                                 dem= c(fun= "quadraticFun", a = 100, b = 5, c = 0),
                                 slope= c(fun= "dnorm", mean = 6000, sd= 3000),
                                 canopy= c(fun= "dnorm", mean = 5, sd= 7))

#create a first species from the parameters with big worldclim data
my_first_species <- generateSpFromFun(raster.stack = my_variables[[c("ai","bio4", "bio6", "bio15","gdd3","pop2000", "pop5000", "rad", "dem", "slope", "canopy")]],
                                      formula = "ai * bio4 + bio6+ bio15 +2*gdd3 +pop2000 + pop5000 + rad + dem + slope + canopy/2 ", 
                                      #bio1 2 times more important than bio2
                                      parameters = my_parameters,
                                      plot= TRUE)
my_first_species
#convert probability distribution to presence absence
gene_linear2 <- convertToPA(my_first_species, PA.method = "probability",
                            prob.method = "linear", a = 0.8, b = 0)


#test the curves for the variables
x <- seq(-3000, 3000, length.out = 100)

density <- logisticFun(x,beta = -100, alpha = 2500000)
density <-quadraticFun(x, a = 100, b = 5, c = 0)
density<- dnorm(x, mean= 0,sd = 1500000)

# Plot the response curve
plot(x, density, type = "l", xlab = "x", ylab = "Density", main= "min Temperature of coldest")




################################################################################
#Plot from a PCA
################################################################################
#create my set of variables for the PCA
my.stack <- my_variables[[c("ai","bio4", "bio6", "bio15","gdd3","pop2000", "pop5000", "rad", "dem", "slope", "canopy")]]
#creat a spatial distribution from my stack with a wide breadth
generalist <- generateSpFromPCA(raster.stack = my.stack, sample.points = TRUE,
                                niche.breadth = "wide")
year_1<- convertToPA(generalist,
                     PA.method = "probability",
                     prob.method = "logistic",
                     beta = "random", alpha = -0.07,
                     plot = TRUE, species.prevalence = 0.80)

#create a spatial distribution with a mean of 4 and sd 0.5 one the first axis and a mean of 1 and sd 0.5 on the second axis
specialist <- generateSpFromPCA(raster.stack = my.stack, sample.points = TRUE, means = c(4, 1), sds = c(0.5, 0.5))

