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

### check model metrics
print(o.models$metrics)
print(k.models$metrics)


#####  Part 10 ::: get variable importance for each sp. ---------------------------------------------------------------------------------------------
# O. koreanus == bg1_10000
print(o.models$contrib[[2]])

# K. koreana == bg1_10000
print(k.models$contrib[[2]])

# export
write.csv(o.models$contrib[[2]], 'data/varimp/O.koreanus_var.imp.csv')
write.csv(k.models$contrib[[2]], 'data/varimp/K.koreana_var_imp.csv')


#####  Part 11 ::: response curves ---------------------------------------------------------------------------------------------
# function to pull out response data

