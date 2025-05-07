##### Calibrate current ENM with CHELSA data only

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(ENMwrap)
library(dplyr)
library(terra)
library(sf)


##### Part 1 ::: prep data ---------------------------------------------------------------------------------------------

## load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv') %>% dplyr::select('long', 'lat')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv') %>% dplyr::select('long', 'lat')

head(o.occs)
head(k.occs)

## load bg
# set1
bg1_5000 <- read.csv('data/bg/CHELSA/set1/bg1_5000.csv') %>% dplyr::select('long', 'lat')
bg1_10000 <- read.csv('data/bg/CHELSA/set1/bg1_10000.csv') %>% dplyr::select('long', 'lat')
bg1_15000 <- read.csv('data/bg/CHELSA/set1/bg1_15000.csv') %>% dplyr::select('long', 'lat')

# set2
bg2_5000 <- read.csv('data/bg/CHELSA/set2/bg2_5000.csv') %>% dplyr::select('long', 'lat')
bg2_10000 <- read.csv('data/bg/CHELSA/set2/bg2_10000.csv') %>% dplyr::select('long', 'lat')
bg2_15000 <- read.csv('data/bg/CHELSA/set2/bg2_15000.csv') %>% dplyr::select('long', 'lat')

## load envs
envs <- rast(list.files(path = 'data/masked/CHELSA', pattern = '.bil', full.names = T))
envs <- terra::subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15'))

print(envs)
plot(envs[[1]])


#####  Part 2 ::: model testing  ---------------------------------------------------------------------------------------------

### test models for O. koreanus
# run
o_chelsa_clim_opt <- test_multibg(taxon.name = 'O.koreanus',
                                  occs = o.occs,
                                  envs = envs,
                                  bg.list = list(bg1_5000, bg1_10000, bg1_15000, 
                                                 bg2_5000, bg2_10000, bg2_15000),
                                  tune.args = list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                                                   rm = seq(0.5,5, by = 0.5)),
                                  partitions = 'checkerboard',
                                  partition.settings = list(aggregation.factor = c(4,4)),
                                  type = 'type1')

# look at optimal param combinations
print(o_chelsa_clim_opt$metrics)

# export metrics
write.csv(o_chelsa_clim_opt$metrics, 'tuning_experiments/preds_addn/O.koreanus/CHELSA_clim_only/metrics/O.koreanus_chelsa_addn_clim_only_metrics.csv')

# look at predictions
print(o_chelsa_clim_opt$preds)
names(o_chelsa_clim_opt$preds) = c('bg1_5000', 'bg1_10000', 'bg1_15000', 'bg2_5000', 'bg2_10000', 'bg2_15000')
plot(o_chelsa_clim_opt$preds)

# export predictions
for (i in 1:nlayers(o_chelsa_clim_opt$preds)) {
  writeRaster(o_chelsa_clim_opt$preds[[i]],
              paste0('tuning_experiments/preds_addn/O.koreanus/CHELSA_clim_only/preds/', names(o_chelsa_clim_opt$preds)[i], '.tif'),
              overwrite = T)
}


### test models for K. koreana
# run
k_chelsa_clim_opt <- test_multibg(taxon.name = 'K.koreana',
                                  occs = k.occs,
                                  envs = envs,
                                  bg.list = list(bg1_5000, bg1_10000, bg1_15000,
                                                 bg2_5000, bg2_10000, bg2_15000),
                                  tune.args = list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                                                   rm = seq(0.5,5, by = 0.5)),
                                  partitions = 'checkerboard',
                                  partition.settings = list(aggregation.factor = c(4,4)),
                                  type = 'type1')

# look at optimal param combinations
print(k_chelsa_clim_opt$metrics)

# export metrics
write.csv(k_chelsa_clim_opt$metrics, 'tuning_experiments/preds_addn/K.koreana/CHELSA_clim_only/metrics/K.koreana_chelsa_addn_clim_only_metrics.csv')

# look at predictions
print(k_chelsa_clim_opt$preds)
names(k_chelsa_clim_opt$preds) = c('bg1_5000', 'bg1_10000', 'bg1_15000', 'bg2_5000', 'bg2_10000', 'bg2_15000')
plot(k_chelsa_clim_opt$preds)

# export predictions
for (i in 1:nlayers(k_chelsa_clim_opt$preds)) {
  writeRaster(k_chelsa_clim_opt$preds[[i]],
              paste0('tuning_experiments/preds_addn/K.koreana/CHELSA_clim_only/preds/', names(k_chelsa_clim_opt$preds)[i], '.tif'),
              overwrite = T)
}


#####  Part 3 ::: model hindcasting  ---------------------------------------------------------------------------------------------

### load hindcast layers
# mPWP
mpwp <- rast(list.files(path = 'data/hindcast_layers/processed_CHELSA/mPWP/', pattern = '.bil$', full.names = T))
print(mpwp)
plot(mpwp[[1]])

# MIS19
mis <- rast(list.files(path = 'data/hindcast_layers/processed_CHELSA/MIS19/', pattern = '.bil$', full.names = T))
print(mis)
plot(mis[[1]])

# LIG
lig <- rast(list.files(path = 'data/hindcast_layers/processed_CHELSA/LIG/', pattern = '.bil$', full.names = T))
print(lig)
plot(lig[[1]])

# LGM
lgm <- rast(list.files(path = 'data/hindcast_layers/processed_CHELSA/LGM/', pattern = '.bil$', full.names = T))
print(lgm)
plot(lgm[[1]])

# MH
mh <- rast(list.files(path = 'data/hindcast_layers/processed_CHELSA/MH/', pattern = '.bil$', full.names = T))
print(mh)
plot(mh[[1]])


### O.koreanus hindcast
o_hind_chelsa <- model_predictr(model = o_chelsa_clim_opt$models[[2]],
                                preds.list = list(mpwp, mis, lig, lgm, mh),
                                pred.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'),
                                method = 'single2multi')

print(o_hind_chelsa)
plot(o_hind_chelsa)

# export predictions
for (i in 1:nlayers(o_hind_chelsa)) {
  writeRaster(o_hind_chelsa[[i]], paste0('tuning_experiments/preds_addn/O.koreanus/CHELSA_clim_only/hindcast/', names(o_hind_chelsa)[i], '.tif'), overwrite = T)
}


### K.koreana hindcast
k_hind_chelsa <- model_predictr(model = k_chelsa_clim_opt$models[[2]],
                                preds.list = list(mpwp, mis, lig, lgm, mh),
                                pred.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'),
                                method = 'single2multi')

print(k_hind_chelsa)
plot(k_hind_chelsa)

# export predictions
for (i in 1:nlayers(k_hind_chelsa)) {
  writeRaster(k_hind_chelsa[[i]], paste0('tuning_experiments/preds_addn/K.koreana/CHELSA_clim_only/hindcast/', names(k_hind_chelsa)[i], '.tif'), overwrite = T)
}
