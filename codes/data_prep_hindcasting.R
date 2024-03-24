#################  hindcasting == use layers acquired from PaleoClim DB
# clear working environment
rm(list = ls(all.names = T))
gc()

# turn off scientific notation
options(scipen = 999)

# load packages
library(raster)
library(dplyr)

##### Part 14 ::: hindcasting data prep  ---------------------------------------------------------------------------------------------
# clipping extent
ext <- c(120, 135, 33, 44)


##### hindcast layer prep
## load current envs == climate only
envs <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil$', full.names = T))
envs <- raster::stack(subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
envs <- raster::setMinMax(envs)
print(envs)

## unit check function
unitCheck <- function(ref.env, proj.env, n) {
  require(dplyr)
  require(dismo)
  
  results <- list()
  
  name.ref.env <- sort(names(ref.env))
  name.proj.env <- sort(names(proj.env))
  
  pts <- randomPoints(mask = ref.env[[1]], n = n) %>% as.data.frame() 
  
  check <- for (i in 1:length(name.ref.env)) {
    ref.ex <- raster::extract(ref.env, pts) %>% na.omit() %>% as.data.frame() 
    proj.ex <- raster::extract(proj.env, pts) %>% na.omit() %>% as.data.frame() 
    
    ref.max <- max(ref.ex[[name.ref.env[[i]]]])
    proj.max <- max(proj.ex[[name.proj.env[[i]]]])
    
    results[[i]] <- data.frame(ref.max, proj.max)
  }
  results <- dplyr::bind_rows(results)
  results$var.name = name.ref.env
  return(results)
}


#### Mid-Pliocene Warm Period (mPWP) == 3.205 Ma
mpwp <- raster::stack(list.files(path = 'data/hindcast_layers/mPWP_v1_r2_5m/2_5min', pattern = '.tif$', full.names = T))
names(mpwp) = gsub('_', '', names(mpwp))

mpwp <- raster::crop(mpwp, ext)
mpwp <- raster::stack(subset(mpwp, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15'))) 
plot(mpwp[[1]])

#### Marine Isotope Stage 19 (MIS19) == 787 Ka
mis <- raster::stack(list.files(path = 'data/hindcast_layers/MIS19_v1_r2_5m/2_5min', pattern = '.tif$', full.names = T))
names(mis) = gsub('_', '', names(mis))

mis <- raster::crop(mis, ext)
mis <- raster::stack(subset(mis, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
plot(mis[[1]])

#### Last Interglacial (LIG) == 130 Ka
lig <- raster::stack(list.files(path = 'data/hindcast_layers/LIG_v1_2_5m/2_5min', pattern = '.tif$', full.names = T))
names(lig) = gsub('_', '', names(lig))

lig <- raster::crop(lig, ext)
lig <- raster::stack(subset(lig, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
plot(lig[[1]])

#### Last Glacial Maximum (LGM) == 21 Ka
lgm <- raster::stack(list.files(path = 'data/hindcast_layers/chelsa_LGM_v1_2B_r2_5m/2_5min', pattern = '.tif$', full.names = T))
names(lgm) = gsub('_', '', names(lgm))

lgm <- raster::crop(lgm, ext)
lgm <- raster::stack(subset(lgm, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
plot(lgm[[1]])

#### Mid-Holocene
mh <- raster::stack(list.files(path = 'data/hindcast_layers/MH_v1_2_5m/2_5min', pattern = '.tif$', full.names = T))
names(mh) = gsub('_', '', names(mh))

mh <- raster::crop(mh, ext)
mh <- raster::stack(subset(mh, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
plot(mh[[1]])


#### check units....
unitCheck(ref.env = envs, proj.env = mpwp, n = 10000)  # mPWP
unitCheck(ref.env = envs, proj.env = mis, n = 10000)  # MIS19
unitCheck(ref.env = envs, proj.env = lig, n = 10000)  # LIG
unitCheck(ref.env = envs, proj.env = lgm, n = 10000)  # LGM
unitCheck(ref.env = envs, proj.env = mh, n = 10000)  # MH

#### ....and "fix" the rasters as needed /// divide temp related variables by 10 == bio1 & bio4
#### NOTE == the max values from "unitCheck" output and printed max value of the full RasterStack can differ
#### this is simply because the "unitCheck" function uses n random points sampled across the reference extent to speed up the calculation! 

# mPWP
mpwp.bio1 <- mpwp[['bio1']]/10
mpwp.bio4 <- mpwp[['bio4']]/10

mpwp <- dropLayer(mpwp, c('bio1', 'bio4'))
mpwp <- raster::stack(mpwp.bio1, mpwp.bio4, mpwp)
print(mpwp)

unitCheck(ref.env = envs, proj.env = mpwp, n = 10000)

# MIS19
mis.bio1 <- mis[['bio1']]/10
mis.bio4 <- mis[['bio4']]/10

mis <- dropLayer(mis, c('bio1', 'bio4'))
mis <- raster::stack(mis.bio1, mis.bio4, mis)
print(mis)

unitCheck(ref.env = envs, proj.env = mis, n = 10000)

# LIG
lig.bio1 <- lig[['bio1']]/10
lig.bio4 <- lig[['bio4']]/10

lig <- dropLayer(lig, c('bio1', 'bio4'))
lig <- raster::stack(lig.bio1, lig.bio4, lig)
print(lig)

unitCheck(ref.env = envs, proj.env = lig, n = 10000)

# LGM
lgm.bio1 <- lgm[['bio1']]/10
lgm.bio4 <- lgm[['bio4']]/10

lgm <- dropLayer(lgm, c('bio1', 'bio4'))
lgm <- raster::stack(lgm.bio1, lgm.bio4, lgm)
print(lgm)

unitCheck(ref.env = envs, proj.env = lgm, n = 10000)

# MH
mh.bio1 <- mh[['bio1']]/10
mh.bio4 <- mh[['bio4']]/10

mh <- dropLayer(mh, c('bio1', 'bio4'))
mh <- raster::stack(mh.bio1, mh.bio4, mh)
print(mh)

unitCheck(ref.env = envs, proj.env = mh, n = 10000)


#### downscale to 1km using bilinear interpolation == fact is 5 since the original resolution (5km) is 5 times greater than the desired resolution
# mPWP
mpwp <- raster::disaggregate(mpwp, fact = 5, method = 'bilinear')
print(mpwp)
plot(mpwp[[1]])

# MIS19
mis <- raster::disaggregate(mis, fact = 5, method = 'bilinear')
print(mis)
plot(mis[[1]])

# LIG
lig <- raster::disaggregate(lig, fact = 5, method = 'bilinear')
print(lig)
plot(lig[[1]])

# LGM
lgm <- raster::disaggregate(lgm, fact = 5, method = 'bilinear')
print(lgm)
plot(lgm[[1]])

# MH
mh <- raster::disaggregate(mh, fact = 5, method = 'bilinear')
print(mh)
plot(mh[[1]])


#### export processed layers
# mPWP
for (i in 1:nlayers(mpwp)) {
  r <- mpwp[[i]]
  name <- paste0('data/hindcast_layers/processed/mPWP/', names(mpwp)[i], '.bil')
  writeRaster(r, filename = name, overwrite = T)
}

# MIS19
for (i in 1:nlayers(mis)) {
  r <- mis[[i]]
  name <- paste0('data/hindcast_layers/processed/MIS19/', names(mis)[i], '.bil')
  writeRaster(r, filename = name, overwrite = T)
}

# LIG
for (i in 1:nlayers(lig)) {
  r <- lig[[i]]
  name <- paste0('data/hindcast_layers/processed/lig/', names(lig)[i], '.bil')
  writeRaster(r, filename = name, overwrite = T)
}

# LGM
for (i in 1:nlayers(lgm)) {
  r <- lgm[[i]]
  name <- paste0('data/hindcast_layers/processed/lgm/', names(lgm)[i], '.bil')
  writeRaster(r, filename = name, overwrite = T)
}

# MH
for (i in 1:nlayers(mh)) {
  r <- mh[[i]]
  name <- paste0('data/hindcast_layers/processed/mh/', names(mh)[i], '.bil')
  writeRaster(r, filename = name, overwrite = T)
}
