##### Run Multivariate Environmental Similarity Surface (MESS) // with WorldClim data ------------------------------------------------------------

# clean working environments
rm(list = ls(all.names = T))
gc()

# load packages
library(raster)
library(dismo)
library(dplyr)

##### part 1 ::: load occs  ------------------------------------------------------------
#o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
#k.occs <- read.csv('data/occs/Karsenia_koreana.csv')

#head(o.occs)
#head(k.occs)


##### part 2 ::: reference env  ------------------------------------------------------------

# reference env rasters
envs <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil$', full.names = T))
envs <- raster::stack(subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
print(envs)
plot(envs[[1]])

# extract reference env values
#o.ref <- raster::extract(envs, o.occs[, c(2,3)]) %>% as.data.frame()
#k.ref <- raster::extract(envs, k.occs[, c(2,3)]) %>% as.data.frame()

#head(o.ref)
#head(k.ref)


##### part 3 ::: load hindcast layers ------------------------------------------------------------

# mid-Pliocene Warm Period (mPWP)
mpwp <- raster::stack(list.files(path = 'data/hindcast_layers/processed_WorldClim/mPWP', pattern = '.bil$', full.names = T))
print(mpwp)
plot(mpwp[[1]])

# Marine Isotope Stage 19 (MIS19)
mis <- raster::stack(list.files(path = 'data/hindcast_layers/processed_WorldClim/MIS19', pattern = '.bil$', full.names = T))
print(mis)
plot(mis[[1]])

# Last Interglacial (LIG)
lig <- raster::stack(list.files(path = 'data/hindcast_layers/processed_WorldClim/LIG', pattern = '.bil$', full.names = T))
print(lig)
plot(lig[[1]])

# Last Glacial Maximum (LGM)
lgm <- raster::stack(list.files(path = 'data/hindcast_layers/processed_WorldClim/LGM', pattern = '.bil$', full.names = T))
print(lgm)
plot(lgm[[1]])

# Mid-Holocene (MH)
mh <- raster::stack(list.files(path = 'data/hindcast_layers/processed_WorldClim/MH', pattern = '.bil$', full.names = T))
print(mh)
plot(mh[[1]])


##### part 4 ::: run MESS ------------------------------------------------------------

### O.koreanus
# MH
#o.mess.mh <- dismo::mess(x = mh, v = o.ref, full = F)
#plot(o.mess.mh$mess)

# LGM
#o.mess.lgm <- dismo::mess(x = lgm, v = o.ref, full = F)
#plot(o.mess.lgm$mess)

# LIG
#o.mess.lig <- dismo::mess(x = lig, v = o.ref, full = F)
#plot(o.mess.lig$mess)

# MIS19
#o.mess.mis <- dismo::mess(x = mis, v = o.ref, full = F)

# mPWP
#o.mess.mpwp <- dismo::mess(x = mpwp, v = o.ref, full = F)


### K.koreana
# MH
#k.mess.mh <- dismo::mess(x = mh, v = k.ref, full = F)
#plot(k.mess.mh$mess)

# LGM
#k.mess.lgm <- dismo::mess(x = lgm, v = k.ref, full = F)

# LIG
#k.mess.lig <- dismo::mess(x = lig, v = k.ref, full = F)

# MIS19
#k.mess.mis <- dismo::mess(x = mis, v = k.ref, full = F)

# mPWP
#k.mess.mpwp <- dismo::mess(x = mpwp, v = k.ref, full = F)


#--------------------------------------------------------------------------------------------------------------------------------

### use ntbox to run mess

# automate
get_mess <- function(m_stack, g_stack_list, names_mess) {
  require(ntbox)
  require(raster)
  
  mess_out <- list()
  
  for (i in 1:length(g_stack_list)) {
    run_mess <- ntb_mess(M_stack = m_stack, G_stack = g_stack_list[[i]])
    mess_out[[i]] <- run_mess
  }
  mess_stack <- raster::stack(mess_out)
  names(mess_stack) = names_mess
  
  return(mess_stack)
}

# run
worldclim.mess <- get_mess(m_stack = envs, g_stack_list = list(mh, lgm, lig, mis, mpwp), 
                           names_mess = c('MH','LGM','LIG','MIS','mPWP'))

# export WorldClim MESS
for (i in 1:nlayers(worldclim.mess)) {
  r <- worldclim.mess[[i]]
  writeRaster(r, filename = paste0('output_other/WorldClim/clim_only/MESS/', names(worldclim.mess)[i], '_ntb_mess.tif'), overwrite = T)
}

