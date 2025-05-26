############   make binary maps from additional predictions

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(ENMwrap)
library(terra)


#####  load occurrence points

# O.koreanus
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
nrow(o.occs)
head(o.occs)

# K.koreana
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')
nrow(k.occs)
head(k.occs)


#####  CHELSA-optimized climate-only  ------------------------------------------

### O.koreanus
# load prediction
o.clim.ch <- rast('tuning_experiments/preds_addn/O.koreanus/CHELSA_clim_only/preds/bg1_10000.tif')
plot(o.clim.ch)

# get threshold
o.clim.ch.th <- get_thresh(preds = raster::raster(o.clim.ch), occs.list = list(o.occs[, c(2,3)]), type = 'p10')
print(o.clim.ch.th)

# get binary
o.clim.ch.bin <- bin_maker(preds = raster::raster(o.clim.ch), th = o.clim.ch.th[[1]])
plot(o.clim.ch.bin)

# export
writeRaster(o.clim.ch.bin, 'tuning_experiments/preds_addn/O.koreanus/bin_both/O.koreanus_addn_clim_only_bin.tif', overwrite = T)


### K.koreana
# load prediction
k.clim.ch <- rast('tuning_experiments/preds_addn/K.koreana/CHELSA_clim_only/preds/bg1_10000.tif')
plot(k.clim.ch)

# get threshold
k.clim.ch.th <- get_thresh(preds = raster::raster(k.clim.ch), occs.list = list(k.occs[, c(2,3)]), type = 'p10')
print(k.clim.ch.th)

# get binary
k.clim.ch.bin <- bin_maker(preds = raster::raster(k.clim.ch), th = k.clim.ch.th[[1]])
plot(k.clim.ch.bin)

# export
writeRaster(k.clim.ch.bin, 'tuning_experiments/preds_addn/K.koreana/bin_both/K.koreana_addn_clim_only_bin.tif', overwrite = T)



#####  CHELSA-optimized full  ------------------------------------------

### O.koreanus
# load prediction
o.full.ch <- rast('tuning_experiments/preds_addn/O.koreanus/CHELSA_full/preds/bg1_10000.tif')
plot(o.full.ch)

# get threshold
o.full.ch.th <- get_thresh(preds = raster::raster(o.full.ch), occs.list = list(o.occs[, c(2,3)]), type = 'p10')
print(o.full.ch.th)

# get binary
o.full.ch.bin <- bin_maker(preds = raster::raster(o.full.ch), th = o.full.ch.th[[1]])
plot(o.full.ch.bin)

# export
writeRaster(o.full.ch.bin, 'tuning_experiments/preds_addn/O.koreanus/bin_both/O.koreanus_addn_full_bin.tif', overwrite = T)


### K.koreana
# load prediction
k.full.ch <- rast('tuning_experiments/preds_addn/K.koreana/CHELSA_full/preds/bg1_10000.tif')
plot(k.full.ch)

# get threshold
k.full.ch.th <- get_thresh(preds = raster::raster(k.full.ch), occs.list = list(k.occs[, c(2,3)]), type = 'p10')
print(k.full.ch.th)

# get binary
k.full.ch.bin <- bin_maker(preds = raster::raster(k.full.ch), th = k.full.ch.th[[1]])
plot(k.full.ch.bin)

# export
writeRaster(k.full.ch.bin, 'tuning_experiments/preds_addn/K.koreana/bin_both/K.koreana_addn_full_bin.tif', overwrite = T)
