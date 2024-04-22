#####  Quantify difference between WorldClim-based and CHELSA-based predictions
# clean up working env
rm(list = ls(all.names = T))
gc()


# load libraries
library(raster)

######### use the formula in Dubos et al. 2023 Biol. Invasions

# function to calculate diff
calc_pred_diff <- function(layer.a, layer.b) {
  sum.a <- cellStats(x = layer.a, stat = 'sum', na.rm = T)
  sum.b <- cellStats(x = layer.b, stat = 'sum', na.rm = T)
  diff <- (abs(sum.a - sum.b)/sum.a) * 100
  
  return(diff)
}

##### O.koreanus
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

# calculate diff
o.cur <- calc_pred_diff(layer.a = o.cur.ch, layer.b = o.cur.wc)
o.mh <- calc_pred_diff(layer.a = o.mh.ch, layer.b = o.mh.wc)
o.lgm <- calc_pred_diff(layer.a = o.lgm.ch, layer.b = o.lgm.wc)
o.lig <- calc_pred_diff(layer.a = o.lig.ch, layer.b = o.lig.wc)
o.mis <- calc_pred_diff(layer.a = o.mis.ch, layer.b = o.mis.wc)
o.mpwp <- calc_pred_diff(layer.a = o.mpwp.ch, layer.b = o.mpwp.wc)


##### K.koreana
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

# calculate diff
k.cur <- calc_pred_diff(layer.a = k.cur.ch, layer.b = k.cur.wc)
k.mh <- calc_pred_diff(layer.a = k.mh.ch, layer.b = k.mh.wc)
k.lgm <- calc_pred_diff(layer.a = k.lgm.ch, layer.b = k.lgm.wc)
k.lig <- calc_pred_diff(layer.a = k.lig.ch, layer.b = k.lig.wc)
k.mis <- calc_pred_diff(layer.a = k.mis.ch, layer.b = k.mis.wc)
k.mpwp <- calc_pred_diff(layer.a = k.mpwp.ch, layer.b = k.mpwp.wc)


