##### Current ENM with WorldClim data 
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


#####  PART 1 ::: load occurrence points ------------------------------------------------------------------------------------------------
# these data are already spatially thinned to 1km 
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')


#####  PART 2 ::: load env data          ------------------------------------------------------------------------------------------------
# load mask polygon 
poly <- rgdal::readOGR('data/polygons/kor_mer.shp')

## climate == WorldClim (1970-2000)
clim <- raster::stack(list.files(path = 'data/WorldClim', pattern = '.tif$', full.names = T))
clim <- raster::crop(clim, extent(poly))
clim <- raster::mask(clim, poly)

plot(clim[[1]])

## topo
topo <- raster::stack(list.files(path = 'data/topo', pattern = '.tif$', full.names = T))
topo <- raster::crop(topo, extent(poly))
topo <- raster::mask(topo, poly)

plot(topo[[1]])

## land cover
land <- raster('data/veg/forest_merged.tif')
land <- raster::crop(land, extent(poly))
land <- raster::mask(land, poly)

plot(land[[1]])

## stack together
envs <- raster::stack(clim, topo, land)
print(envs)


#####  PART 3 ::: background data          ------------------------------------------------------------------------------------------------
### set 1 ::: target group background points == all available amphibian points from the Korean Peninsula

# list of spp 
spplist <- read.csv('data/target_group/baseline/korean_amphibian_splist.csv') %>% as.vector()
print(spplist)

# collect occs
megaSDM::OccurrenceCollection(spplist = spplist$Species,
                              output = 'data/target_group',
                              trainingarea = extent(poly))






