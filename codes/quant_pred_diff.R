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
  diff <- abs(sum.a - sum.b)/sum.a * 100
  
  return(diff)
}

### calculate diff
# O.koreanus
o.cur <- calc_pred_diff(layer.a = o.cur.wc, layer.b = o.cur.ch)
o.mh <- calc_pred_diff(layer.a = o.mh.wc, layer.b = o.mh.ch)
o.lgm <- calc_pred_diff(layer.a = o.lgm.wc, layer.b = o.lgm.ch)
o.lig <- calc_pred_diff(layer.a = o.lig.wc, layer.b = o.lig.ch)
o.mis <- calc_pred_diff(layer.a = o.mis.wc, layer.b = o.mis.ch)
o.mpwp <- calc_pred_diff(layer.a = o.mpwp.wc, layer.b = o.mpwp.ch)

print(round(o.cur, digits = 0))
print(round(o.mh, digits = 0))
print(round(o.lgm, digits = 0))
print(round(o.lig, digits = 0))
print(round(o.mis, digits = 0))
print(round(o.mpwp, digits = 0))

# K.koreana
k.cur <- calc_pred_diff(layer.a = k.cur.wc, layer.b = k.cur.ch)
k.mh <- calc_pred_diff(layer.a = k.mh.wc, layer.b = k.mh.ch)
k.lgm <- calc_pred_diff(layer.a = k.lgm.wc, layer.b = k.lgm.ch)
k.lig <- calc_pred_diff(layer.a = k.lig.wc, layer.b = k.lig.ch)
k.mis <- calc_pred_diff(layer.a = k.mis.wc, layer.b = k.mis.ch)
k.mpwp <- calc_pred_diff(layer.a = k.mpwp.wc, layer.b = k.mpwp.ch)

print(round(k.cur, digits = 0))
print(round(k.mh, digits = 0))
print(round(k.lgm, digits = 0))
print(round(k.lig, digits = 0))
print(round(k.mis, digits = 0))
print(round(k.mpwp, digits = 0))


##### Part 3 ::: use ENMTools ----------------------------------------------------------------------------------------------------

# function
schoener_diff <- function(layer.a, layer.b, verbose = T) {
  sim <- ENMTools::raster.overlap(x = terra::rast(layer.a), y = terra::rast(layer.b), verbose = T)
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

print(round(o.cur2, digits = 0))
print(round(o.mh2, digits = 0))
print(round(o.lgm2, digits = 0))
print(round(o.lig2, digits = 0))
print(round(o.mis2, digits = 0))
print(round(o.mpwp2, digits = 0))

# K.koreana
k.cur2 <- schoener_diff(layer.a = k.cur.ch, layer.b = k.cur.wc, verbose = T)
k.mh2 <- schoener_diff(layer.a = k.mh.ch, layer.b = k.mh.wc, verbose = T)
k.lgm2 <- schoener_diff(layer.a = k.lgm.ch, layer.b = k.lgm.wc, verbose = T)
k.lig2 <- schoener_diff(layer.a = k.lig.ch, layer.b = k.lig.wc, verbose = T)
k.mis2 <- schoener_diff(layer.a = k.mis.ch, layer.b = k.mis.wc, verbose = T)
k.mpwp2 <- schoener_diff(layer.a = k.mpwp.ch, layer.b = k.mpwp.wc, verbose = T)

print(round(k.cur2, digits = 0))
print(round(k.mh2, digits = 0))
print(round(k.lgm2, digits = 0))
print(round(k.lig2, digits = 0))
print(round(k.mis2, digits = 0))
print(round(k.mpwp2, digits = 0))


##### Part 4 ::: export results ----------------------------------------------------------------------------------------------------
# O.koreanus overall
o.diff <- data.frame(round(o.cur, digits = 0), round(o.mh, digits = 0), round(o.lgm, digits = 0), 
                     round(o.lig, digits = 0), round(o.mis, digits = 0), round(o.mpwp, digits = 0))

o.diff$species = 'O.koreanus'
o.diff$method = 'overall'

colnames(o.diff) = c('current', 'MH', 'LGM', 'LIG', 'MIS19', 'mPWP', 'species', 'method')
print(o.diff)

# K.koreana overall
k.diff <- data.frame(round(k.cur, digits = 0), round(k.mh, digits = 0), round(k.lgm, digits = 0), 
                     round(k.lig, digits = 0), round(k.mis, digits = 0), round(k.mpwp, digits = 0))

k.diff$species = 'K.koreana'
k.diff$method = 'overall'

colnames(k.diff) = c('current', 'MH', 'LGM', 'LIG', 'MIS19', 'mPWP', 'species', 'method')
print(k.diff)

# O.koreanus Schoener
o.diff2 <- data.frame(round(o.cur2, digits = 0), round(o.mh2, digits = 0), round(o.lgm2, digits = 0), 
                      round(o.lig2, digits = 0), round(o.mis2, digits = 0), round(o.mpwp2, digits = 0))

o.diff2$species = 'O.koreanus'
o.diff2$method = 'schoener'

colnames(o.diff2) = c('current', 'MH', 'LGM', 'LIG', 'MIS19', 'mPWP', 'species', 'method')
print(o.diff2)

# K.koreana Schoener
k.diff2 <- data.frame(round(k.cur2, digits = 0), round(k.mh2, digits = 0), round(k.lgm2, digits = 0), 
                      round(k.lig2, digits = 0), round(k.mis2, digits = 0), round(k.mpwp2, digits = 0))

k.diff2$species = 'K.koreana'
k.diff2$method = 'schoener'

colnames(k.diff2) = c('current', 'MH', 'LGM', 'LIG', 'MIS19', 'mPWP', 'species', 'method')
print(k.diff2)

# bind the data & export
diff.data <- rbind(o.diff, o.diff2, k.diff, k.diff2)
diff.data <- diff.data[, c(7,8,1,2,3,4,5,6)]
print(diff.data)

write.csv(diff.data, 'clim_source_compare/quant_pred_diff.csv')
