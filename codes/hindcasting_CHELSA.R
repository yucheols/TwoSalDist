#######  hindcasting from the climate-only ENMs
# clean up working env
rm(list = ls(all.names = T))
gc()

### load packages
library(raster)
library(dplyr)
library(ENMeval)
library(ggplot2)


##### Part 16 ::: hindcasting run --------------------------------------------------------------------------------------------------------------------

### load hindcast layers
# mid-Pliocene Warm Period (mPWP)
mpwp <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/mPWP', pattern = '.bil$', full.names = T))
print(mpwp)
plot(mpwp[[1]])

# Marine Isotope Stage 19 (MIS19)
mis <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/MIS19', pattern = '.bil$', full.names = T))
print(mis)
plot(mis[[1]])

# Last Interglacial (LIG)
lig <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/LIG', pattern = '.bil$', full.names = T))
print(lig)
plot(lig[[1]])

# Last Glacial Maximum (LGM)
lgm <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/LGM', pattern = '.bil$', full.names = T))
print(lgm)
plot(lgm[[1]])

# Mid-Holocene (MH)
mh <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/MH', pattern = '.bil$', full.names = T))
print(mh)
plot(mh[[1]])


### load climate-only models
# O. koreanus 
o.models_clim <- readRDS('tuning_experiments/output_model_rds/O_koreanus_clim_only_CHELSA.rds')
print(o.models_clim)

k.models_clim <- readRDS('tuning_experiments/output_model_rds/K_koreana_clim_only_CHELSA.rds')
print(k.models_clim)


##### model prediction to past climatic conditions
# make function for prediction
model_predictr <- function(model, preds.list, pred.names) {
  require(dismo)
  
  output <- list()
  
  for (i in 1:length(preds.list)) {
    make.pred <- dismo::predict(object = model, x = preds.list[[i]], progress = 'text')
    output[[i]] <- make.pred
  }
  output.pred <- raster::stack(output)
  names(output.pred) = pred.names
  return(output.pred)
}

### O.koreanus
# run
o.hinds <- model_predictr(model = o.models_clim$models[[1]], 
                          preds.list = list(mpwp, mis, lig, lgm, mh), 
                          pred.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'))

# print & plot
print(o.hinds)
plot(o.hinds)


### K.koreana
# run
k.hinds <- model_predictr(model = k.models_clim$models[[1]],
                          preds.list = list(mpwp, mis, lig, lgm, mh),
                          pred.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'))

# print & plot
print(k.hinds)
plot(k.hinds)
