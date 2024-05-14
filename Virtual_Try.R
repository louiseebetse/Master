library(virtualspecies)
library(terra)
library(geodata)
library(git2r)
library(sf)

setwd("/Users/louise/Desktop/UNIL/Master/Rstudio/Virtual species")
par(mfrow=c(1,1))

################################################################################
#variables avec bioclim
#worldclim
worldclim<- worldclim_global(var= "bio", res=10, path= "worldclim")
worldclim[[1:4]]
names(worldclim)<-paste0("bio", 1:19)
#plot the raster
plot(worldclim[[c("wc2.1_10m_bio_1", "wc2.1_10m_bio_12")]])

#create the parameters
my_parameters <- formatFunctions(bio1 =c(fun= "dnorm", mean = 25, sd= 5),
                                 bio12 = c(fun= "dnorm", mean = 4000, 
                                           sd= 2000))
#create a first species from the parameters with big worldclim data
my_first_species <- generateSpFromFun(raster.stack = worldclim[[c("bio1", "bio12")]],
                                      formula = "2 * bio1 + bio12", 
                                      #bio1 2 times more important than bio2
                                      parameters = my_parameters,
                                      plot= TRUE)
my_first_species
################################################################################

##variables with swiss bioclim
#Extract data from worldclim for switzerland with the lowest resolution(10)
swissclim<- worldclim_country(country= "Switzerland" ,var= "bio", res=10, path= "worldclim")
swiss_shape<-st_read("swiss_shapefile.gpkg")
swissclim<- mask(swissclim,swiss_shape, inverse= FALSE)
#change the names of the rasters
names(swissclim)<-paste0("bio", 1:19)
plot(swissclim[[c("bio1", "bio6", "bio12","bio14")]])

#parameters for a generalist
gene_para<- formatFunctions(bio1 =c(fun= "dnorm", mean = 20, sd= 6),
                            bio12 = c(fun= "dnorm", mean = 4000, 
                                      sd= 2000),
                            bio6 = c( fun= "logisticFun", beta = 0, alpha = -4),
                            bio14=c(fun= "dnorm", mean = 110, 
                                    sd= 25))
#create the swiss virtual species from the parameters
gene_species <- generateSpFromFun(raster.stack = swissclim[[c("bio1", "bio6", "bio12","bio14")]],
                                      formula = "bio1 + bio12 + 2 * bio6 + bio14", 
                                      #bio6 (coldest t°) 2 times more important than bio2
                                      parameters = gene_para,
                                      plot= TRUE)
#create a distribution with a threshold
gene_threshold <- convertToPA(gene_species, PA.method = "threshold", beta = 0.65)
#create a distribution with a linear method
gene_linear <- convertToPA(gene_species, PA.method = "probability",
                        prob.method = "linear", a = 1, b = 0)
gene_linear2 <- convertToPA(gene_species, PA.method = "probability",
                           prob.method = "linear", a = 0.5, b = 0)
#create a distribution with a logistic method
gene_probability <- convertToPA(gene_species,
                                  PA.method = "probability",
                                  prob.method = "logistic",
                                  beta = 0.65, alpha = -0.07,
                                  plot = TRUE)

############################
#Evolution of the prevalence through years
year_1<- convertToPA(gene_species,
                     PA.method = "probability",
                     prob.method = "logistic",
                     beta = "random", alpha = -0.07,
                     plot = TRUE, species.prevalence = 0.80)

year_2<- convertToPA(gene_species,
                     PA.method = "probability",
                     prob.method = "logistic",
                     beta = "random", alpha = -0.07,
                     plot = TRUE, species.prevalence = 0.70)

year_3<- convertToPA(gene_species,
                     PA.method = "probability",
                     prob.method = "logistic",
                     beta = "random", alpha = -0.07,
                     plot = TRUE, species.prevalence = 0.60)

year_4<- convertToPA(gene_species,
                     PA.method = "probability",
                     prob.method = "logistic",
                     beta = "random", alpha = -0.07,
                     plot = TRUE, species.prevalence = 0.50)
year_1$PA.conversion[5]
year_1$PA.conversion["species.prevalence"]
round(as.numeric(year_1$PA.conversion[5]), 2)


par(mfrow=c(2,2))

plot(year_1$pa.raster, main= paste("Presence of generalist with a prevalence of", round(as.numeric(year_1$PA.conversion[5]), 2)), cex.main= 0.8)
plot(year_2$pa.raster, main= paste("Presence of generalist with a prevalence of", round(as.numeric(year_2$PA.conversion[5]), 2)), cex.main= 0.8)
plot(year_3$pa.raster, main= paste("Presence of generalist with a prevalence of", round(as.numeric(year_3$PA.conversion[5]), 2)), cex.main= 0.8)
plot(year_4$pa.raster, main= paste("Presence of generalist with a prevalence of", round(as.numeric(year_4$PA.conversion[5]), 2)), cex.main= 0.8)


#plot of the distribution of species
plotResponse(gene_species)



#to test the response curves
x <- seq(-20, 50, length.out = 100)

# Calculate the density values using dnorm
density <- logisticFun(x,beta = 0, alpha = -4)

# Plot the response curve
plot(x, density, type = "l", xlab = "x", ylab = "Density", main= "min Temperature of coldest")


##to generate a specialist from a PCA

my.stack <- swissclim[[c("bio2", "bio5", "bio6", "bio12", "bio13", "bio14")]]
specialist <- generateSpFromPCA(raster.stack = my.stack, sample.points = TRUE,
                            niche.breadth = "wide")
specialist <- generateSpFromPCA(raster.stack = my.stack, sample.points = TRUE, means = c(4, 1), sds = c(0.5, 0.5))

#---------------------------------------

#Linear function: formatFunctions(bio1 = c(fun = 'linearFun', a = 1, b = 0))
#Quadratic function: formatFunctions(bio1 = c(fun = 'quadraticFun', a = -1, b = 2, c = 0))
#Bell-shaped function: formatFunctions(bio1 = c(fun = 'dnorm', mean = 25, sd = 5))
#Logistic function: formatFunctions(bio1 = c(fun = 'logisticFun', beta = 150, alpha = -5))
#Normal function defined by extremes: formatFunctions(bio1 = c(fun = 'custnorm', mean = 25, diff = 5, prob = 0.99))
#Beta response function: formatFunctions(bio1 = c(fun = 'betaFun', p1 = 0, p2 = 25, alpha = 0.9, gamma = 0.08))
#other : https://borisleroy.com/virtualspecies_tutorial/02-response.html#how-to-create-and-use-your-own-response-functions

#---------------------------------------
