##### Calibrate current ENM with CHELSA data + topo + vegetation data 

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
bg1_10000 <- read.csv('data/bg/CHELSA/set1/bg1_10000.csv') %>% dplyr::select('long', 'lat')

## load envs
envs <- rast(list.files(path = 'data/masked/CHELSA', pattern = '.bil', full.names = T))
envs <- terra::subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15', 'forest', 'slope'))

print(envs)
plot(envs[[1]])


#####  Part 2 ::: model testing  ---------------------------------------------------------------------------------------------

### test models for O. koreanus
# run
o_chelsa_full_opt <- test_multibg(taxon.name = 'O.koreanus',
                                  occs = o.occs,
                                  envs = envs,
                                  bg.list = list(bg1_10000),
                                  tune.args = list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                                                   rm = seq(0.5,5, by = 0.5)),
                                  partitions = 'checkerboard',
                                  partition.settings = list(aggregation.factor = c(4,4)),
                                  type = 'type1')

# look at optimal param combinations
print(o_chelsa_full_opt$metrics)

# export metrics
write.csv(o_chelsa_full_opt$metrics, 'tuning_experiments/preds_addn/O.koreanus/CHELSA_full/metrics/O.koreanus_chelsa_addn_full_metrics.csv')

# look at predictions
print(o_chelsa_full_opt$preds)
plot(o_chelsa_full_opt$preds)

# export predictions
writeRaster(o_chelsa_full_opt$preds, 'tuning_experiments/preds_addn/O.koreanus/CHELSA_full/preds/bg1_10000.tif', overwrite = T)


### test models for K. koreana
k_chelsa_full_opt <- test_multibg(taxon.name = 'K.koreana',
                                  occs = k.occs,
                                  envs = envs,
                                  bg.list = list(bg1_10000),
                                  tune.args = list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                                                   rm = seq(0.5,5, by = 0.5)),
                                  partitions = 'checkerboard',
                                  partition.settings = list(aggregation.factor = c(4,4)),
                                  type = 'type1')


# look at optimal param combinations
print(k_chelsa_full_opt$metrics)

# export metrics
write.csv(k_chelsa_full_opt$metrics, 'tuning_experiments/preds_addn/K.koreana/CHELSA_full/metrics/K.koreana_chelsa_addn_full_metrics.csv')

# look at predictions
print(k_chelsa_full_opt$preds)
plot(k_chelsa_full_opt$preds)

# export predictions
writeRaster(k_chelsa_full_opt$preds, 'tuning_experiments/preds_addn/K.koreana/CHELSA_full/preds/bg1_10000.tif', overwrite = T)


