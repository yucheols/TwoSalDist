#####  Quantify difference between WorldClim-based and CHELSA-based predictions
# clean up working env
rm(list = ls(all.names = T))
gc()

# load libraries
library(raster)
library(ENMTools)

##### Part 1 ::: load predictions ----------------------------------------------------------------------------------------------------
### O.koreanus
# import WorldClim-based predictions
o.cur.wc <- raster('tuning_experiments/preds/O.koreanus/WorldClim/cont_clim_only/bg1_10000.tif')
o.mh.wc <- raster('hindcast/WorldClim/O_koreanus/O.koreanus_MH.tif')
o.lgm.wc <- raster('hindcast/WorldClim/O_koreanus/O.koreanus_LGM.tif')
o.lig.wc <- raster('hindcast/WorldClim/O_koreanus/O.koreanus_LIG.tif')
o.mis.wc <- raster('hindcast/WorldClim/O_koreanus/O.koreanus_MIS19.tif')
o.mpwp.wc <- raster('hindcast/WorldClim/O_koreanus/O.koreanus_mPWP.tif')

# import CHELSA-based predictions
o.cur.ch <- raster('tuning_experiments/preds/O.koreanus/CHELSA/cont_clim_only/bg1_10000.tif')
o.mh.ch <- raster('hindcast/CHELSA/O_koreanus/O.koreanus_MH.tif')
o.lgm.ch <- raster('hindcast/CHELSA/O_koreanus/O.koreanus_LGM.tif')
o.lig.ch <- raster('hindcast/CHELSA/O_koreanus/O.koreanus_LIG.tif')
o.mis.ch <- raster('hindcast/CHELSA/O_koreanus/O.koreanus_MIS19.tif')
o.mpwp.ch <- raster('hindcast/CHELSA/O_koreanus/O.koreanus_mPWP.tif')


### K.koreana
# import WorldClim-based predictions
k.cur.wc <- raster('tuning_experiments/preds/K.koreana/WorldClim/cont_clim_only/bg1_10000.tif')
k.mh.wc <- raster('hindcast/WorldClim/K_koreana/K.koreana_MH.tif')
k.lgm.wc <- raster('hindcast/WorldClim/K_koreana/K.koreana_LGM.tif')
k.lig.wc <- raster('hindcast/WorldClim/K_koreana/K.koreana_LIG.tif')
k.mis.wc <- raster('hindcast/WorldClim/K_koreana/K.koreana_MIS19.tif')
k.mpwp.wc <- raster('hindcast/WorldClim/K_koreana/K.koreana_mPWP.tif')

# import CHELSA-based predictions
k.cur.ch <- raster('tuning_experiments/preds/K.koreana/CHELSA/cont_clim_only/bg1_10000.tif')
k.mh.ch <- raster('hindcast/CHELSA/K_koreana/K.koreana_MH.tif')
k.lgm.ch <- raster('hindcast/CHELSA/K_koreana/K.koreana_LGM.tif')
k.lig.ch <- raster('hindcast/CHELSA/K_koreana/K.koreana_LIG.tif')
k.mis.ch <- raster('hindcast/CHELSA/K_koreana/K.koreana_MIS19.tif')
k.mpwp.ch <- raster('hindcast/CHELSA/K_koreana/K.koreana_mPWP.tif')


##### Part 2 ::: use the formula in Dubos et al. 2023 Biol. Invasions -------------------------------------------------------------------------------------------------

# function to calculate diff
calc_pred_diff <- function(layer.a, layer.b) {
  sum.a <- cellStats(x = layer.a, stat = 'sum', na.rm = T)
  sum.b <- cellStats(x = layer.b, stat = 'sum', na.rm = T)
  diff <- (abs(sum.a - sum.b)/sum.a) * 100
  
  return(diff)
}

### calculate diff
# O.koreanus
o.cur <- calc_pred_diff(layer.a = o.cur.ch, layer.b = o.cur.wc)
o.mh <- calc_pred_diff(layer.a = o.mh.ch, layer.b = o.mh.wc)
o.lgm <- calc_pred_diff(layer.a = o.lgm.ch, layer.b = o.lgm.wc)
o.lig <- calc_pred_diff(layer.a = o.lig.ch, layer.b = o.lig.wc)
o.mis <- calc_pred_diff(layer.a = o.mis.ch, layer.b = o.mis.wc)
o.mpwp <- calc_pred_diff(layer.a = o.mpwp.ch, layer.b = o.mpwp.wc)

print(o.cur)
print(o.mh)
print(o.lgm)
print(o.lig)
print(o.mis)
print(o.mpwp)

# K.koreana
k.cur <- calc_pred_diff(layer.a = k.cur.ch, layer.b = k.cur.wc)
k.mh <- calc_pred_diff(layer.a = k.mh.ch, layer.b = k.mh.wc)
k.lgm <- calc_pred_diff(layer.a = k.lgm.ch, layer.b = k.lgm.wc)
k.lig <- calc_pred_diff(layer.a = k.lig.ch, layer.b = k.lig.wc)
k.mis <- calc_pred_diff(layer.a = k.mis.ch, layer.b = k.mis.wc)
k.mpwp <- calc_pred_diff(layer.a = k.mpwp.ch, layer.b = k.mpwp.wc)

print(k.cur)
print(k.mh)
print(k.lgm)
print(k.lig)
print(k.mis)
print(k.mpwp)


##### Part 3 ::: use ENMTools ----------------------------------------------------------------------------------------------------

# function
schoener_diff <- function(layer.a, layer.b, verbose = T) {
  sim <- ENMTools::raster.overlap(x = layer.a, y = layer.b, verbose = T)
  sim_diff <- (1 - sim$D) * 100
  
  return(sim_diff)
}

# O.koreanus
o.cur2 <- schoener_diff(layer.a = o.cur.ch, layer.b = o.cur.wc, verbose = T)
o.mh2 <- schoener_diff(layer.a = o.mh.ch, layer.b = o.mh.wc, verbose = T)
o.lgm2 <- schoener_diff(layer.a = o.lgm.ch, layer.b = o.lgm.wc, verbose = T)
o.lig2 <- schoener_diff(layer.a = o.lig.ch, layer.b = o.lig.wc, verbose = T)
o.mis2 <- schoener_diff(layer.a = o.mis.ch, layer.b = o.mis.wc, verbose = T)
o.mpwp2 <- schoener_diff(layer.a = o.mpwp.ch, layer.b = o.mpwp.wc, verbose = T)

print(o.cur2)
print(o.mh2)
print(o.lgm2)
print(o.lig2)
print(o.mis2)
print(o.mpwp2)

# K.koreana
k.cur2 <- schoener_diff(layer.a = k.cur.ch, layer.b = k.cur.wc, verbose = T)
k.mh2 <- schoener_diff(layer.a = k.mh.ch, layer.b = k.mh.wc, verbose = T)
k.lgm2 <- schoener_diff(layer.a = k.lgm.ch, layer.b = k.lgm.wc, verbose = T)
k.lig2 <- schoener_diff(layer.a = k.lig.ch, layer.b = k.lig.wc, verbose = T)
k.mis2 <- schoener_diff(layer.a = k.mis.ch, layer.b = k.mis.wc, verbose = T)
k.mpwp2 <- schoener_diff(layer.a = k.mpwp.ch, layer.b = k.mpwp.wc, verbose = T)

print(k.cur2)
print(k.mh2)
print(k.lgm2)
print(k.lig2)
print(k.mis2)
print(k.mpwp2)


##### Part 4 ::: export results ----------------------------------------------------------------------------------------------------
quant_pred_diff <- data.frame()




