#####  post modeling stuff == model eval, var contributions, resp curves, etc.
#####  this continues from the "model_WorldClim.R" workflow
getwd()

# clean up working env
rm(list = ls(all.names = T))
gc()

# load packages
library(raster)
library(dplyr)
library(ENMeval)


#####  Part 10 ::: get variable importance for each sp. ---------------------------------------------------------------------------------------------
# O. koreanus
print(o.models$contrib[[2]])

# K. koreana
print(k.models$contrib[[1]])
