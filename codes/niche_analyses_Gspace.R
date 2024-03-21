##### Niche analyses in G-space using ENMTools
# setwd
setwd('niche_analyses/G-space')

# clean up working env
rm(list = ls(all.names = T))
gc()

# load packages
library(ENMTools)
library(raster)
library(ggplot2)
library(dplyr)


##### Step 1 ::: prep data ----------------------------------------------------------------------------------------------------

### occurrence data
occs.sp1 <- read.csv('occs/Onychodactylus_koreanus.csv')
occs.sp2 <- read.csv('occs/Karsenia_koreana.csv')

head(occs.sp1)
head(occs.sp2)

### background data 
bg.sp1 <- read.csv('bg/set1/bg1_10000.csv')
bg.sp2 <- read.csv('bg/set1/bg1_10000.csv')

### envs data
envs <- raster::stack(list.files(path = 'masked/WorldClim', pattern = '.bil$', full.names = T))
envs <- raster::stack(subset(envs, c('bio1', 'bio3', 'bio4', 'bio12', 'bio14', 'forest', 'slope')))

print(envs)
plot(envs[[1]])

### species ranges == 50km radius dissolved circular buffer around each occ point 
sp1.range <- rgdal::readOGR('range_poly/O.koreanus/o_range.shp')
sp2.range <- rgdal::readOGR('range_poly/K.koreana/k_range.shp')


##### Step 2 ::: create enmtools.species objects -------------------------------------------------------------------------------
# sp1
sp1 <- enmtools.species(range = raster::mask(envs[[1]], sp1.range),
                        presence.points = occs.sp1[, c(2,3)],
                        background.points = bg.sp1[, c(2,3)],
                        species.name = 'O.koreanus')

# sp2
sp2 <- enmtools.species(range = raster::mask(envs[[1]], sp2.range),
                        presence.points = occs.sp2[, c(2,3)],
                        background.points = bg.sp2[, c(2,3)],
                        species.name = 'K.koreana')


##### Step 3 ::: identity test -------------------------------------------------------------------------------
id.test <- identity.test(species.1 = sp1, species.2 = sp2, env = envs, type = 'mx', nreps = 100, 
                         bg.source = 'points', low.memory = T, rep.dir = 'id_test', verbose = T, clamp = T)


##### Step 4 ::: background test -------------------------------------------------------------------------------
# O.koreanus vs K.koreana background
bg.test1 <- background.test(species.1 = sp1, species.2 = sp2, env = envs, type = 'mx', nreps = 100, test.type = 'asymmetric',
                            bg.source = 'points', low.memory = T, rep.dir = 'bg_test1', verbose = T, clamp = T)

# K.koreana vs O.koreanus background
bg.test2 <- background.test(species.1 = sp2, species.2 = sp1, env = envs, type = 'mx', nreps = 100, test.type = 'asymmetric',
                            bg.source = 'points', low.memory = T, rep.dir = 'bg_test2', verbose = T, clamp = T)


## Save outputs
#saveRDS(id.test, 'output/id_test.rds')
#saveRDS(bg.test1, 'output/bg_test1.rds')
#saveRDS(bg.test2, 'output/bg_test2.rds')


##### Step 5 ::: plot results -------------------------------------------------------------------------------
