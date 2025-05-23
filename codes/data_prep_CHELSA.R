##### Data prep for current ENM with CHELSA data 
# clean up working env
rm(list = ls(all.names = T))
gc()

# set seed 
set.seed(111)

# prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

# load packages
library(raster)
library(plyr)
library(dplyr) 
library(megaSDM)
library(readr)
library(blockCV)
library(sf)
library(tmap)

#####  PART 1 ::: load occurrence points ------------------------------------------------------------------------------------------------
# these data are already spatially thinned to 1km 
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')


#####  PART 2 ::: load env data          ------------------------------------------------------------------------------------------------
# load mask polygon 
poly <- rgdal::readOGR('data/polygons/kor_mer.shp')

## climate == CHELSA (1979-2013)
clim <- raster::stack(list.files(path = 'data/CHELSA/30sec', pattern = '.tif$', full.names = T))
clim <- raster::crop(clim, extent(poly))
clim <- raster::mask(clim, poly)

plot(clim[[1]])

names(clim) = gsub('_', '', names(clim))
print(clim)

# drop layers not represented in paleoclim layers
clim <- dropLayer(clim, c('bio2', 'bio3','bio5', 'bio6', 'bio7'))  

## topo
topo <- raster::stack(list.files(path = 'data/topo', pattern = '.tif$', full.names = T))
topo <- raster::crop(topo, extent(poly))
topo <- raster::mask(topo, poly)

plot(topo[[1]])

## land cover
land <- raster('data/veg/forest.tif')
land <- raster::crop(land, extent(poly))
land <- raster::mask(land, poly)

plot(land[[1]])

## stack together
envs <- raster::stack(clim, topo, land)
print(envs)

## export masked layers == .bil format
for (i in 1:nlayers(envs)) {
  layer <- envs[[i]]
  name <- paste0('data/masked/CHELSA/', names(envs)[i], '.bil')
  writeRaster(layer, filename = name, overwrite = T)
}

## create import shortcut
# CAUTION ::: For CHELSA layers.....importing the .bil files through the shortcut below will throw errors at the background sampling step for whatever reason
# CAUTION ::: For CHELSA data.....use the original .tif import /// use .bil files for layer selection

#envs <- raster::stack(list.files(path = 'data/masked/CHELSA', pattern = '.bil$', full.names = T))
#plot(envs[[1]])


#####  PART 3 ::: background data          ------------------------------------------------------------------------------------------------
### use the bg data sampled in "data_prep_WorldClim" code

# import thinned target group points
#targ.pts <- read.csv('data/target_group/thinned/targ.pts.csv') %>% dplyr::select(2,3)
#head(targ.pts)

###  make density raster for Set 1
#targ.ras1 <- rasterize(targ.pts, envs, 1)
#plot(targ.ras1)

#targ.pres1 <- which(values(targ.ras1) == 1)
#targ.pres.locs1 <- coordinates(targ.ras1)[targ.pres1, ]

#targ.dens1 <- MASS::kde2d(targ.pres.locs1[,1], targ.pres.locs1[,2],
#                          n = c(nrow(targ.ras1), ncol(targ.ras1)),
#                          lims = c(extent(envs)[1], extent(envs)[2], extent(envs)[3], extent(envs)[4]))

#targ.dens.ras1 <- raster(targ.dens1, envs)
#targ.dens.ras1 <- resample(targ.dens.ras1, envs)

#bias.layer1 <- raster::mask(targ.dens.ras1, poly)
#plot(bias.layer1)


### sample background points at three different sample sizes == 5,000 // 10,000 // 15,000
# n = 5000
#bg1_5000 <- xyFromCell(bias.layer1,
#                       sample(which(!is.na(values(subset(envs, 1)))), 5000,
#                              prob = values(bias.layer1)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

#colnames(bg1_5000) = c('long', 'lat')
#points(bg1_5000, col = 'blue')

#write.csv(bg1_5000, 'data/bg/CHELSA/set1/bg1_5000.csv')


# n = 10000
#bg1_10000 <- xyFromCell(bias.layer1,
#                        sample(which(!is.na(values(subset(envs, 1)))), 10000,
#                               prob = values(bias.layer1)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

#colnames(bg1_10000) = c('long', 'lat')
#points(bg1_10000, col = 'blue')

#write.csv(bg1_10000, 'data/bg/CHELSA/set1/bg1_10000.csv')


# n = 15000
#bg1_15000 <- xyFromCell(bias.layer1,
#                        sample(which(!is.na(values(subset(envs, 1)))), 15000,
#                              prob = values(bias.layer1)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

#colnames(bg1_15000) = c('long', 'lat')
#points(bg1_15000, col = 'blue')

#write.csv(bg1_15000, 'data/bg/CHELSA/set1/bg1_15000.csv')


### set 2 ::: pooled occurrence points
#targ.ras2 <- rasterize(rbind(o.occs[, c(2,3)], k.occs[, c(2,3)]), envs, 1)
#plot(targ.ras2)

#targ.pres2 <- which(values(targ.ras2) == 1)
#targ.pres.locs2 <- coordinates(targ.ras2)[targ.pres2, ]

#targ.dens2 <- MASS::kde2d(targ.pres.locs2[,1], targ.pres.locs2[,2],
#                          n = c(nrow(targ.ras2), ncol(targ.ras2)),
#                          lims = c(extent(envs)[1], extent(envs)[2], extent(envs)[3], extent(envs)[4]))

#targ.dens.ras2 <- raster(targ.dens2, envs)
#targ.dens.ras2 <- resample(targ.dens.ras2, envs)

#bias.layer2 <- raster::mask(targ.dens.ras2, poly)
#plot(bias.layer2)

### sample background points at three different sample sizes == 5,000 // 10,000 // 15,000
# n = 5000
#bg2_5000 <- xyFromCell(bias.layer2,
#                       sample(which(!is.na(values(subset(envs, 1)))), 5000,
#                              prob = values(bias.layer2)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

#colnames(bg2_5000) = c('long', 'lat')
#points(bg2_5000, col = 'green')

#write.csv(bg2_5000, 'data/bg/CHELSA/set2/bg2_5000.csv')


# n = 10000
#bg2_10000 <- xyFromCell(bias.layer2,
#                        sample(which(!is.na(values(subset(envs, 1)))), 10000,
#                               prob = values(bias.layer2)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

#colnames(bg2_10000) = c('long', 'lat')
#points(bg2_10000, col = 'green')

#write.csv(bg2_10000, 'data/bg/CHELSA/set2/bg2_10000.csv')


# n = 15000
#bg2_15000 <- xyFromCell(bias.layer2,
#                        sample(which(!is.na(values(subset(envs, 1)))), 15000,
#                               prob = values(bias.layer2)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

#colnames(bg2_15000) = c('long', 'lat')
#points(bg2_15000, col = 'green')

#write.csv(bg2_15000, 'data/bg/CHELSA/set2/bg2_15000.csv')


### use fixed bg
bg1_10000 <- read.csv('data/bg/bg1_10000.csv')

#####  PART 4 ::: select environmental data    ------------------------------------------------------------------------------------------------
# use ntbox
ntbox::run_ntbox()

### Pearson |r| > 0.8 removed ==  bio1 bio4 bio12 bio13 bio14 bio15 forest slope  
envs <- raster::stack(subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15', 'forest', 'slope')))

#print(envs)
plot(envs[[1]])


#####  PART 5 ::: make cross validation folds...if you want to specify them yourself  ----------------------------------------------------------------------
# make "user-specified" folds

### Step 1 === get block size
# automate block size calc
block_size <- function(points, raster, num_sample, crs) {
  
  block_size_list <- list()
  
  block_prep <- for (i in 1:length(points)) {
    points[[i]]$occ <- rep(0, nrow(points[[i]]))
    pts.sf <- sf::st_as_sf(points[[i]], coords = c('long', 'lat'), crs = crs)
    auto  <- blockCV::cv_spatial_autocor(x = pts.sf, r = raster, num_sample = num_sample, column = 'occ', progress = T)
    
    block_size_list[[i]] <- auto$range
  }
  return(block_size_list)
}

# get block size
get.block.size <- block_size(points = list(bg1_5000, bg1_10000, bg1_15000, bg2_5000, bg2_10000, bg2_15000),
                             raster = envs[[1:6]], num_sample = 5000, crs = 4326)

get.block.size[1:3] # bg set 1
get.block.size[4:6] # bg set 2


### Step 2 === create folds
# automate fold generation
fold_maker <- function(occs, bg.list, envs, k, block.size) {
  
  # output
  occ.fold.out <- list()
  bg.fold.out <- list()
  
  # make loop
  fold.maker <- for (i in 1:length(bg.list)) {
    
    # bind occ & bg ::: both datasets should have matching column names
    occs_bg <- rbind(occs, bg.list[[i]])
    occs_bg$occ <- c(rep(1, nrow(occs)), rep(0, nrow(bg.list[[i]])))
    
    # extract envs values ::: check NA & remove them
    occs_bg_env <- raster::extract(envs, occs_bg[, c(1,2)]) %>% as.data.frame()
    occs_bg_env <- cbind(occs_bg_env, occs_bg) %>% na.omit()
    
    # select point columns
    #pts_cols <- occs_bg_env[, c(10:12)]  # deprecated
    pts_cols <- occs_bg_env[, c('long', 'lat', 'occ')]
    
    # convert to simple features
    pts_sf <- sf::st_as_sf(pts_cols, coords = c('long', 'lat'), crs = 4326)
    
    # create spatial folds
    spat_folds <- blockCV::cv_spatial(x = pts_sf, column = 'occ', r = envs, k = k, hexagon = T, flat_top = F,
                                      size = block.size[[i]], selection = 'random', iteration = 100, seed = 333, 
                                      progress = T, plot = F)
    
    # segregate into occ & bg folds
    folds_ids <- as.data.frame(spat_folds$folds_ids)
    colnames(folds_ids) = 'fold_id'
    
    folds <- cbind(pts_cols, folds_ids)
    
    occ_folds <- folds %>% dplyr::filter(occ == 1)
    bg_folds <- folds %>% dplyr::filter(occ == 0)
    
    occ.fold.out[[i]] <- occ_folds
    bg.fold.out[[i]] <- bg_folds
    
  }
  
  return(list(occ.fold.out, bg.fold.out))
  
}

######## make folds per species
#### O. koreanus
### O. koreanus -- bg set1 -- 4 fold
o.folds.bg1 <- fold_maker(occs = o.occs[, c(2,3)], 
                          bg.list = list(bg1_5000[, c('long', 'lat')], bg1_10000[, c('long', 'lat')], bg1_15000[, c('long', 'lat')]), 
                          envs = envs[[1:6]], k = 4, block.size = get.block.size[1:3])

# occ folds
o.occ.folds.bg1 <- o.folds.bg1[[1]]

# bg folds
o.bg.folds.bg1 <- o.folds.bg1[[2]]


### O. koreanus -- bg set2 -- 4 fold
o.folds.bg2 <- fold_maker(occs = o.occs[, c(2,3)],
                          bg.list = list(bg2_5000[, c('long', 'lat')], bg2_10000[, c('long', 'lat')], bg2_15000[, c('long', 'lat')]),
                          envs = envs[[1:6]], k = 4, block.size = get.block.size[4:6])

# occ folds
o.occ.folds.bg2 <- o.folds.bg2[[1]]

# bg folds
o.bg.folds.bg2 <- o.folds.bg2[[2]]



#### K. koreana
### K. koreana -- bg set1 -- 2 fold
k.folds.bg1 <- fold_maker(occs = k.occs[, c(2,3)], 
                          bg.list = list(bg1_5000[, c('long', 'lat')], bg1_10000[, c('long', 'lat')], bg1_15000[, c('long', 'lat')]), 
                          envs = envs[[1:6]], k = 2, block.size = get.block.size[1:3])

# occ folds
k.occ.folds.bg1 <- k.folds.bg1[[1]]

# bg folds
k.bg.folds.bg1 <- k.folds.bg1[[2]]


### O. koreanus -- bg set2 -- 4 fold
k.folds.bg2 <- fold_maker(occs = k.occs[, c(2,3)],
                          bg.list = list(bg2_5000[, c('long', 'lat')], bg2_10000[, c('long', 'lat')], bg2_15000[, c('long', 'lat')]),
                          envs = envs[[1:6]], k = 2, block.size = get.block.size[4:6])

# occ folds
k.occ.folds.bg2 <- k.folds.bg2[[1]]

# bg folds
k.bg.folds.bg2 <- k.folds.bg2[[2]]



### Step 3 === compile folds
# O. koreanus == 4-fold
o.folds <- list(list(occs.grp = o.occ.folds.bg1[[1]]$fold_id, bg.grp = o.bg.folds.bg1[[1]]$fold_id),
                list(occs.grp = o.occ.folds.bg1[[2]]$fold_id, bg.grp = o.bg.folds.bg1[[2]]$fold_id),
                list(occs.grp = o.occ.folds.bg1[[3]]$fold_id, bg.grp = o.bg.folds.bg1[[3]]$fold_id),
                list(occs.grp = o.occ.folds.bg2[[1]]$fold_id, bg.grp = o.bg.folds.bg2[[1]]$fold_id),
                list(occs.grp = o.occ.folds.bg2[[2]]$fold_id, bg.grp = o.bg.folds.bg2[[2]]$fold_id),
                list(occs.grp = o.occ.folds.bg2[[3]]$fold_id, bg.grp = o.bg.folds.bg2[[3]]$fold_id))

print(o.folds)


# K. koreana == 2-fold
k.folds <- list(list(occs.grp = k.occ.folds.bg1[[1]]$fold_id, bg.grp = k.bg.folds.bg1[[1]]$fold_id),
                list(occs.grp = k.occ.folds.bg1[[2]]$fold_id, bg.grp = k.bg.folds.bg1[[2]]$fold_id),
                list(occs.grp = k.occ.folds.bg1[[3]]$fold_id, bg.grp = k.bg.folds.bg1[[3]]$fold_id),
                list(occs.grp = k.occ.folds.bg2[[1]]$fold_id, bg.grp = k.bg.folds.bg2[[1]]$fold_id),
                list(occs.grp = k.occ.folds.bg2[[2]]$fold_id, bg.grp = k.bg.folds.bg2[[2]]$fold_id),
                list(occs.grp = k.occ.folds.bg2[[3]]$fold_id, bg.grp = k.bg.folds.bg2[[3]]$fold_id))

print(k.folds)


### Step 4 === check fold assignments == the number and assignment of folds should be EXACTLY THE SAME between the occ & bg
#             otherwise ENMeval will throw the following error == task 1 failed - "cannot evaluate a model without absence and presence data that are not NA

### O. koreanus
# occs
for (i in 1:length(o.folds)) {
  folds <- unique(sort(o.folds[[i]]$occs.grp))
  print(folds)
}

# bg
for (i in 1:length(o.folds)) {
  folds <- unique(sort(o.folds[[i]]$bg.grp))
  print(folds)
}


### K. koreana
for (i in 1:length(k.folds)) {
  folds <- unique(sort(k.folds[[i]]$occs.grp))
  print(folds)
}

# bg
for (i in 1:length(k.folds)) {
  folds <- unique(sort(k.folds[[i]]$bg.grp))
  print(folds)
}


### save folds
#saveRDS(o.folds, 'data/folds/CHELSA/O.koreanus_user_folds.rds')
#saveRDS(k.folds, 'data/folds/CHELSA/K.koreana_user_folds.rds')



#####  PART 5b ::: explore data partitioning methods  ----------------------------------------------------------------------
#### data extraction for downstream use

## O.koreanus
o.occs.z <- cbind(o.occs[, -1], raster::extract(envs, o.occs[, -1]))
o.bg.z <- cbind(bg1_10000[, -1], raster::extract(envs, bg1_10000[, -1]))

head(o.occs.z)
head(o.bg.z)

## K.koreana
k.occs.z <- cbind(k.occs[, -1], raster::extract(envs, k.occs[, -1]))
k.bg.z <- cbind(bg1_10000[, -1], raster::extract(envs, bg1_10000[, -1]))

head(k.occs.z)
head(k.occs.z)

########  try checkerboard2 

#### O.koreanus
o.cb2 <- ENMeval::get.checkerboard2(occs = o.occs[, -1], bg = bg1_10000[, -1], envs = envs[[1:6]], aggregation.factor = c(4,4))
ENMeval::evalplot.grps(pts = o.occs[, -1], pts.grp = o.cb2$occs.grp, envs = envs)
ENMeval::evalplot.grps(pts = bg1_10000[, -1], pts.grp = o.cb2$bg.grp, envs = envs)

# check points per bin
table(o.cb2$occs.grp)
table(o.cb2$bg.grp)

# check environmental similarity between bins
ENMeval::evalplot.envSim.hist(sim.type = 'mess', ref.data = 'occs', occs.z = o.occs.z, bg.z = o.bg.z, occs.grp = o.cb2$occs.grp, bg.grp = o.cb2$bg.grp) # occs
ENMeval::evalplot.envSim.hist(sim.type = 'mess', ref.data = 'bg', occs.z = o.occs.z, bg.z = o.bg.z, occs.grp = o.cb2$occs.grp, bg.grp = o.cb2$bg.grp)   # bg


#### K.koreana
k.cb2 <- ENMeval::get.checkerboard2(occs = k.occs[, -1], bg = bg1_10000[, -1], envs = envs[[1:6]], aggregation.factor = c(4,4))
ENMeval::evalplot.grps(pts = k.occs[, -1], pts.grp = k.cb2$occs.grp, envs = envs)
ENMeval::evalplot.grps(pts = bg1_10000[, -1], pts.grp = k.cb2$bg.grp, envs = envs)

# check points per bin
table(k.cb2$occs.grp)
table(k.cb2$bg.grp)

# check environmental similarity between bins
ENMeval::evalplot.envSim.hist(sim.type = 'mess', ref.data = 'occs', occs.z = k.occs.z, bg.z = k.bg.z, occs.grp = k.cb2$occs.grp, bg.grp = k.cb2$bg.grp) # occs
ENMeval::evalplot.envSim.hist(sim.type = 'mess', ref.data = 'bg', occs.z = k.occs.z, bg.z = k.bg.z, occs.grp = k.cb2$occs.grp, bg.grp = k.cb2$bg.grp)   # bg
