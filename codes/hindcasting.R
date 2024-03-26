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
mpwp <- raster::stack(list.files(path = 'data/hindcast_layers/processed/mPWP', pattern = '.bil$', full.names = T))
print(mpwp)
plot(mpwp[[1]])

# Marine Isotope Stage 19 (MIS19)
mis <- raster::stack(list.files(path = 'data/hindcast_layers/processed/MIS19', pattern = '.bil$', full.names = T))
print(mis)
plot(mis[[1]])

# Last Interglacial (LIG)
lig <- raster::stack(list.files(path = 'data/hindcast_layers/processed/LIG', pattern = '.bil$', full.names = T))
print(lig)
plot(lig[[1]])

# Last Glacial Maximum (LGM)
lgm <- raster::stack(list.files(path = 'data/hindcast_layers/processed/LGM', pattern = '.bil$', full.names = T))
print(lgm)
plot(lgm[[1]])

# Mid-Holocene (MH)
mh <- raster::stack(list.files(path = 'data/hindcast_layers/processed/MH', pattern = '.bil$', full.names = T))
print(mh)
plot(mh[[1]])


### load climate-only models
# O. koreanus 
o.models_clim <- readRDS('tuning_experiments/output_model_rds/O_koreanus_clim_only_WorldClim.rds')
print(o.models_clim)

k.models_clim <- readRDS('tuning_experiments/output_model_rds/K_koreana_clim_only_WorldClim.rds')
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


### export hindcast rasters
# O.koreanus
for (i in 1:nlayers(o.hinds)) {
  r <- o.hinds[[i]]
  name <- paste0('hindcast/O_koreanus/O.koreanus_', names(o.hinds)[i], '.tif')
  writeRaster(r, filename = name, overwrite = T)
}

# K.koreana
for (i in 1:nlayers(k.hinds)) {
  r <- k.hinds[[i]]
  name <- paste0('hindcast/K_koreana/K.koreana_', names(k.hinds)[i], '.tif')
  writeRaster(r, filename = name, overwrite = T)
}


##### Part 17 ::: MESS --------------------------------------------------------------------------------------------------------------------


##### Part 18 ::: hindcast binary for LGM and MH == prep for dispersal analyses using SDMtoolbox in ArcGIS --------------------------------

### O.koreanus == clim only p10 == 0.3719397
# LGM
o.lgm.bin <- ecospat::ecospat.binary.model(Pred = terra::rast(o.hinds$LGM), Threshold = 0.3719397) %>% raster()
plot(o.lgm.bin)

# MH
o.mh.bin <- ecospat::ecospat.binary.model(Pred = terra::rast(o.hinds$MH), Threshold = 0.3719397) %>% raster()
plot(o.mh.bin)

# export
writeRaster(o.lgm.bin, 'dispersal_corridors/rasters/O.koreanus/O.koreanus_LGM_bin.tif', overwrite = T)
writeRaster(o.mh.bin, 'dispersal_corridors/rasters/O.koreanus/O.koreanus_MH_bin.tif', overwrite = T)


### K.koreana == clim only p10 == 0.4518484
# LGM
k.lgm.bin <- ecospat::ecospat.binary.model(Pred = terra::rast(k.hinds$LGM), Threshold = 0.4518484) %>% raster()
plot(k.lgm.bin)

# MH
k.mh.bin <- ecospat::ecospat.binary.model(Pred = terra::rast(k.hinds$MH), Threshold = 0.4518484) %>% raster()
plot(k.mh.bin)

# export
writeRaster(k.lgm.bin, 'dispersal_corridors/rasters/K.koreana/K.koreana_LGM_bin.tif', overwrite = T)
writeRaster(k.mh.bin, 'dispersal_corridors/rasters/K.koreana/K.koreana_MH_bin.tif', overwrite = T)


##### Part 19 :::  compare env values between time ----------------------------------------------------------------------




##### Part 20 :::  optional == visualize climate trends through time // from mPWP to current --------------------------------------------------------------
##### annual mean temp and annual precip

## function to automate data formatting
env.val.get <- function(env.list, var, time, type) {
  require(dplyr)
  
  output <- list()
  
  if (type == 'max'){
    for (i in 1:length(env.list)) {
      get.df <- env.list[[i]][[var]] %>% as.data.frame() %>% na.omit()
      get.df.max <- max(get.df[[var]])
      df.max <- data.frame(val = get.df.max, time = time[[i]])
      
      output[[i]] <- df.max
      output[[i]]$type = type
      output[[i]]$var = var
    }
  }
  else if (type == 'mean') {
    for (i in 1:length(env.list)) {
      get.df <- env.list[[i]][[var]] %>% as.data.frame() %>% na.omit()
      get.df.mean <- mean(get.df[[var]])
      df.mean <- data.frame(val = get.df.mean, time = time[[i]])
      
      output[[i]] <- df.mean
      output[[i]]$type = type
      output[[i]]$var = var
    }
  }
  else if (type == 'min') {
    for (i in 1:length(env.list)) {
      get.df <- env.list[[i]][[var]] %>% as.data.frame() %>% na.omit()
      get.df.min <- min(get.df[[var]])
      df.min <- data.frame(val = get.df.min, time = time[[i]])
      
      output[[i]] <- df.min
      output[[i]]$type = type
      output[[i]]$var = var
    }
  }
  output.all <- dplyr::bind_rows(output)
  return(output.all)
}


### load current values
cur <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil$', full.names = T))
cur <- raster::stack(subset(cur, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))  
print(cur)
plot(cur[[1]])

### bio1 ::: Annual Mean Temp
# min
min.vals.tmp <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, cur), var = 'bio1', 
                            time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'min')

print(min.vals.tmp)

# mean
mean.vals.tmp <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, cur), var = 'bio1', 
                             time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'mean')

print(mean.vals.tmp)

# max
max.vals.tmp <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, cur), var = 'bio1', 
                            time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'max')

print(max.vals.tmp)


### bio13 ::: 
# min
min.vals.pr <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, cur), var = 'bio13', 
                           time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'min')

print(min.vals.pr)

# mean
mean.vals.pr <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, cur), var = 'bio13', 
                            time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'mean')

print(mean.vals.pr)

# max
max.vals.pr <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, cur), var = 'bio13', 
                           time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'max')

print(max.vals.pr)


### bind
env.vals <- rbind(min.vals.tmp, mean.vals.tmp, max.vals.tmp, min.vals.pr, mean.vals.pr, max.vals.pr)
print(env.vals)

### reorder variable levels 
env.vals$time = factor(env.vals$time, levels = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'))
env.vals$type = factor(env.vals$type, levels = c('min', 'mean', 'max'))

### rename variables
env.vals$var <- dplyr::recode_factor(env.vals$var, 'bio1' = 'Bio1 (°C)', 'bio13' = 'Bio13 (mm)')
env.vals$type <- dplyr::recode_factor(env.vals$type, 'min' = 'Min', 'mean' = 'Mean', 'max' = 'Max')

### plot
env.vals %>%
  ggplot(aes(x = time, y = val, group = type, color = type)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ var, scales = 'free') +
  xlab('Time period') + ylab('Values') +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = 'top')
