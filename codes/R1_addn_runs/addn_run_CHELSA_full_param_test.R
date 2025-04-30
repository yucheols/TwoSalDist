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
## load mask polygon for later plotting
poly <- st_read('data/polygons/kor_mer.shp')

## load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv') %>% select('long', 'lat')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv') %>% select('long', 'lat')

head(o.occs)
head(k.occs)

## load bg
# set1
bg1_5000 <- read.csv('data/bg/CHELSA/set1/bg1_5000.csv') %>% select('long', 'lat')
bg1_10000 <- read.csv('data/bg/CHELSA/set1/bg1_10000.csv') %>% select('long', 'lat')
bg1_15000 <- read.csv('data/bg/CHELSA/set1/bg1_15000.csv') %>% select('long', 'lat')

# set2
bg2_5000 <- read.csv('data/bg/CHELSA/set2/bg2_5000.csv') %>% select('long', 'lat')
bg2_10000 <- read.csv('data/bg/CHELSA/set2/bg2_10000.csv') %>% select('long', 'lat')
bg2_15000 <- read.csv('data/bg/CHELSA/set2/bg2_15000.csv') %>% select('long', 'lat')

## load envs
envs <- rast(list.files(path = 'data/masked/CHELSA', pattern = '.bil', full.names = T))
envs <- terra::subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15', 'forest', 'slope'))

print(envs)
plot(envs[[1]])


#####  Part 15 ::: model testing  ---------------------------------------------------------------------------------------------

### test models for O. koreanus
o_chelsa_full_opt <- test_multibg(taxon.name = 'O.koreanus',
                                  occs = o.occs,
                                  envs = envs,
                                  bg.list = list(bg1_5000, bg1_10000, bg1_15000,
                                                 bg2_5000, bg2_10000, bg2_15000),
                                  tune.args = list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                                                   rm = seq(0.5,5, by = 0.5)),
                                  partitions = 'checkerboard',
                                  partition.settings = list(aggregation.factor = c(4,4)),
                                  type = 'type1')



### test models for K. koreana
k_chelsa_full_opt <- test_multibg(taxon.name = 'K.koreana',
                                  occs = k.occs,
                                  envs = envs,
                                  bg.list = list(bg1_5000, bg1_10000, bg1_15000,
                                                 bg2_5000, bg2_10000, bg2_15000),
                                  tune.args = list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                                                   rm = seq(0.5,5, by = 0.5)),
                                  partitions = 'checkerboard',
                                  partition.settings = list(aggregation.factor = c(4,4)),
                                  type = 'type1')