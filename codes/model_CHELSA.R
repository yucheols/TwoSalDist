##### Current ENM with CHELSA data 
# clean up working env
rm(list = ls(all.names = T))

# set seed 
set.seed(321)

# prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

# load packages
library(raster)

#####  PART 1 ::: load occurrence points ------------------------------------------------------------------------------------------------
# these data are already spatially thinned to 1km 
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')

#####  PART 2 ::: load env data          ------------------------------------------------------------------------------------------------
# load mask polygon 
poly <- rgdal::readOGR('data/polygons/kor_mer.shp')

## climate == CHELSA (1979-2013)
clim <- raster::stack(list.files(path = 'data/CHELSA', pattern = '.tif$', full.names = T))
clim <- raster::mask(clim, poly)

plot(clim[[1]])



