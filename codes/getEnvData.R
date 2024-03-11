#####  get PaleoClim data & CHELSA data using the package rpaleoclim //// WorldClim data are downloaded directly from the website
library(rpaleoclim)
library(terra)

#####  get CHELSA current data == 1979 - 2013
ext <- c(124.1824, 130.9404, 33.11208, 43.00605)

# get data
pc_data <- paleoclim(period = c('cur'), resolution = '30s', region = ext, as = 'raster', cache_path = 'data/CHELSA')
