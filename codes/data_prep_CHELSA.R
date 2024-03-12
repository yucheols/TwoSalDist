##### Data prep for current ENM with CHELSA data 
# clean up working env
rm(list = ls(all.names = T))
gc()

# set seed 
set.seed(321)

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
library(usdm)


#####  PART 1 ::: load occurrence points ------------------------------------------------------------------------------------------------
# these data are already spatially thinned to 1km 
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')


#####  PART 2 ::: load env data          ------------------------------------------------------------------------------------------------
# load mask polygon 
poly <- rgdal::readOGR('data/polygons/kor_mer.shp')

## climate == CHELSA (1979-2013)
clim <- raster::stack(list.files(path = 'data/CHELSA', pattern = '.tif$', full.names = T))
clim <- raster::crop(clim, extent(poly))
clim <- raster::mask(clim, poly)

plot(clim[[1]])

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
envs <- raster::stack(list.files(path = 'data/masked/CHELSA', pattern = '.bil$', full.names = T))
plot(envs[[1]])


#####  PART 3 ::: background data          ------------------------------------------------------------------------------------------------
### set 1 ::: target group background points == all available amphibian points from the Korean Peninsula
# list of spp 
spplist <- read.csv('data/target_group/baseline/korean_amphibian_splist.csv') %>% as.vector()
print(spplist)

# collect target group species points
megaSDM::OccurrenceCollection(spplist = spplist$Species,
                              output = 'data/target_group',
                              trainingarea = extent(poly))

# compile points
targ.pts <- list.files(path = 'data/target_group', pattern = '.csv', full.names = T) %>%
  lapply(read_csv) %>%
  rbind.fill %>%
  as.data.frame() %>%
  dplyr::select(4,6,5)

colnames(targ.pts) = c('species', 'long', 'lat')
head(targ.pts)


# add on occurrence points based on 4th NES amphibian datatpoints and Borzee et al. 2021
nes_nk <- read.csv('data/target_group/baseline/target_group_NES_NK_SK.csv') 
head(nes_nk)

# combine the two
targ.pts <- rbind(targ.pts, nes_nk)

# thin
targ.pts <- SDMtune::thinData(coords = targ.pts[, c(2,3)], env = terra::rast(envs[[1]]), x = 'long', y = 'lat', verbose = T, progress = T)

###  make density raster for Set 1
targ.ras1 <- rasterize(targ.pts, envs, 1)
plot(targ.ras1)

targ.pres1 <- which(values(targ.ras1) == 1)
targ.pres.locs1 <- coordinates(targ.ras1)[targ.pres1, ]

targ.dens1 <- MASS::kde2d(targ.pres.locs1[,1], targ.pres.locs1[,2],
                          n = c(nrow(targ.ras1), ncol(targ.ras1)),
                          lims = c(extent(envs)[1], extent(envs)[2], extent(envs)[3], extent(envs)[4]))

targ.dens.ras1 <- raster(targ.dens1, envs)
targ.dens.ras1 <- resample(targ.dens.ras1, envs)

bias.layer1 <- raster::mask(targ.dens.ras1, poly)
plot(bias.layer1)


### sample background points at three different sample sizes == 5,000 // 10,000 // 15,000
# n = 5000
bg1_5000 <- xyFromCell(bias.layer1,
                       sample(which(!is.na(values(subset(envs, 1)))), 5000,
                              prob = values(bias.layer1)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

colnames(bg1_5000) = c('long', 'lat')
points(bg1_5000, col = 'blue')

write.csv(bg1_5000, 'data/bg/set1/bg1_5000.csv')


# n = 10000
bg1_10000 <- xyFromCell(bias.layer1,
                        sample(which(!is.na(values(subset(envs, 1)))), 10000,
                               prob = values(bias.layer1)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

colnames(bg1_10000) = c('long', 'lat')
points(bg1_10000, col = 'blue')

write.csv(bg1_10000, 'data/bg/set1/bg1_10000.csv')


# n = 15000
bg1_15000 <- xyFromCell(bias.layer1,
                        sample(which(!is.na(values(subset(envs, 1)))), 15000,
                               prob = values(bias.layer1)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

colnames(bg1_15000) = c('long', 'lat')
points(bg1_15000, col = 'blue')

write.csv(bg1_15000, 'data/bg/set1/bg1_15000.csv')


### set 2 ::: pooled occurrence points
targ.ras2 <- rasterize(rbind(o.occs[, c(2,3)], k.occs[, c(2,3)]), envs, 1)
plot(targ.ras2)

targ.pres2 <- which(values(targ.ras2) == 1)
targ.pres.locs2 <- coordinates(targ.ras2)[targ.pres2, ]

targ.dens2 <- MASS::kde2d(targ.pres.locs2[,1], targ.pres.locs2[,2],
                          n = c(nrow(targ.ras2), ncol(targ.ras2)),
                          lims = c(extent(envs)[1], extent(envs)[2], extent(envs)[3], extent(envs)[4]))

targ.dens.ras2 <- raster(targ.dens2, envs)
targ.dens.ras2 <- resample(targ.dens.ras2, envs)

bias.layer2 <- raster::mask(targ.dens.ras2, poly)
plot(bias.layer2)

### sample background points at three different sample sizes == 5,000 // 10,000 // 15,000
# n = 5000
bg2_5000 <- xyFromCell(bias.layer2,
                       sample(which(!is.na(values(subset(envs, 1)))), 5000,
                              prob = values(bias.layer2)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

colnames(bg2_5000) = c('long', 'lat')
points(bg2_5000, col = 'green')

write.csv(bg2_5000, 'data/bg/set2/bg2_5000.csv')


# n = 10000
bg2_10000 <- xyFromCell(bias.layer2,
                        sample(which(!is.na(values(subset(envs, 1)))), 10000,
                               prob = values(bias.layer2)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

colnames(bg2_10000) = c('long', 'lat')
points(bg2_10000, col = 'green')

write.csv(bg2_10000, 'data/bg/set2/bg2_10000.csv')


# n = 15000
bg2_15000 <- xyFromCell(bias.layer2,
                        sample(which(!is.na(values(subset(envs, 1)))), 15000,
                               prob = values(bias.layer2)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame()

colnames(bg2_15000) = c('long', 'lat')
points(bg2_15000, col = 'green')

write.csv(bg2_15000, 'data/bg/set2/bg2_15000.csv')


#####  PART 4 ::: select environmental data    ------------------------------------------------------------------------------------------------
# extract 50000 random points across the extent
pts <- dismo::randomPoints(mask = envs[[1]], n = 50000) %>% as.data.frame()

# extract raster values
vals <- raster::extract(envs, pts)

# vifstep
vifstep(vals, th = 10, maxobservations = 50000)


#####  PART 5 ::: make cross validation folds      ---------------------------------------------------------------------------------------------