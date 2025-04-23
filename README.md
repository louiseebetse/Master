# Master

There are three important programs :
  one_full_go
  One_full_go_random
  Convert_to_matrix
  WS_rarity

## one_full_go
Allows from real presences, environmental rasters to create suitability map, populations and observations in the three different trends (increase, decrease, stable)
Then put them at a 5km resolution and save them again
It also create the table with the value of the populations and observations at each time step

## Convert to matrix
from the observations it is going to create a sp_matrix, a sampling effort matrix, a suitability matrix and the occupancy matrix from the Infoflora model. In it it uses the function occupancy_function2 this is the corrected occupancy function.
HOWEVER it does not take into account the WS map
At the end of this file there is a function to create the graph with all the populations trends and all the infoflora models.

## One_full_go_random
Allows from environmental rasters to create random virtual species. Thus there is a random suitability map, populations and observations in the three different trends (increase, decrease, stable)
Then put them at a 5km resolution and save them again
It also create the table with the value of the populations and observations at each time step
from the observations it is going to create a sp_matrix, a sampling effort matrix, a suitability matrix and the occupancy matrix from the Infoflora model. In it it uses the function occupancy_function2 this is the corrected occupancy function.
HOWEVER it does not take into account the WS map
At the end of this file there is a function to create the graph with all the populations trends and all the infoflora models.

## WS_rarity
From the virtual populations this script creates a map with the WS partitioning of switzerland, and apply if there is no observation, few observation or a lot of observations in every WS shape.
