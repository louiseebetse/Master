library(virtualspecies)
library(terra)
library(geodata)
library(git2r)
library(sf)
library(gtools) # V.3.9.4
setwd("/Users/louise/Desktop/UNIL/Master/Rstudio/Virtual species")

my_variables <- mixedsort(list.files(path = "variables",
                                     pattern = NULL,
                                     full.names = TRUE))

my_variables<- rast(my_variables)
cat(crs(my_variables))
names(my_variables) = c("ai","bio4", "bio6", "bio15","gdd3","pop2000", "pop5000", "rad", "dem", "slope", "canopy")
plot(my_variables)
my_variables

#create the parameters
my_parameters <- formatFunctions(ai = c(fun= "dnorm", mean = 25, sd= 5),
                                 bio4 =c(fun= "dnorm", mean = 25, sd= 5),
                                 bio6 =c(fun= "dnorm", mean = 25, sd= 5),
                                 bio15 = c(fun= "dnorm", mean = 4000, sd= 2000),
                                 gdd3 =c(fun= "dnorm", mean = 25, sd= 5),
                                 pop2000= c(fun= "dnorm", mean = 25, sd= 5),
                                 pop5000= c(fun= "dnorm", mean = 25, sd= 5),
                                 rad= c(fun= "dnorm", mean = 25, sd= 5),
                                 dem= c(fun= "dnorm", mean = 25, sd= 5),
                                 slope= c(fun= "dnorm", mean = 25, sd= 5),
                                 canopy= c(fun= "dnorm", mean = 25, sd= 5))
#create a first species from the parameters with big worldclim data
my_first_species <- generateSpFromFun(raster.stack = worldclim[[c("bio1", "bio12")]],
                                      formula = "2 * bio1 + bio12", 
                                      #bio1 2 times more important than bio2
                                      parameters = my_parameters,
                                      plot= TRUE)
my_first_species


x <- seq(-20, 50, length.out = 100)

# Calculate the density values using dnorm
density <- logisticFun(x,beta = 0, alpha = -4)

# Plot the response curve
plot(x, density, type = "l", xlab = "x", ylab = "Density", main= "min Temperature of coldest")
