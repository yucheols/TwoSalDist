##### Calibrate current ENM with WorldClim data only

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(raster)
library(dplyr)
library(ENMeval)
library(rasterVis)
library(ggplot2)

##### Part 6 ::: prep data ---------------------------------------------------------------------------------------------
## load mask polygon for later plotting
poly <- rgdal::readOGR('data/polygons/kor_mer.shp')

## load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv') %>% dplyr::select('long', 'lat')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv') %>% dplyr::select('long', 'lat')

head(o.occs)
head(k.occs)

## load bg
# set1
bg1_5000 <- read.csv('data/bg/set1/bg1_5000.csv') %>% dplyr::select('long', 'lat')
bg1_10000 <- read.csv('data/bg/set1/bg1_10000.csv') %>% dplyr::select('long', 'lat')
bg1_15000 <- read.csv('data/bg/set1/bg1_15000.csv') %>% dplyr::select('long', 'lat')

# set2
bg2_5000 <- read.csv('data/bg/set2/bg2_5000.csv') %>% dplyr::select('long', 'lat')
bg2_10000 <- read.csv('data/bg/set2/bg2_10000.csv') %>% dplyr::select('long', 'lat')
bg2_15000 <- read.csv('data/bg/set2/bg2_15000.csv') %>% dplyr::select('long', 'lat')

## load envs
envs <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil', full.names = T))
envs <- raster::stack(subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
print(envs)
plot(envs[[1]])


##### Part 7 ::: Model tuning ---------------------------------------------------------------------------------------------
# automate model tuning 
# type 1 == minimum or.10p.avg as primary criterion // type 2 == delta.AICc <= 2 as primary criterion 
test_models <- function(taxon.name, occs, envs, bg.list, tune.args, partitions, partition.settings = NULL, user.grp = NULL, type) {
  output <- list()
  models <- list()
  preds <- list()
  contrib <- list()
  
  if (type == 'type1') {
    for (i in 1:length(bg.list)) {
      
      # make models
      eval <- ENMeval::ENMevaluate(taxon.name = taxon.name, occs = occs, envs = envs, bg = bg.list[[i]], 
                                   tune.args = tune.args, partitions = partitions, partition.settings = partition.settings, 
                                   user.grp = user.grp[[i]], doClamp = T, algorithm = 'maxent.jar', parallel = T, parallelType = 'doSNOW')
      
      # get results
      eval.res <- ENMeval::eval.results(eval)
      
      # get optimal parameter combinations
      opt.param <- eval.res %>% dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
        dplyr::filter(auc.diff.avg == min(auc.diff.avg)) %>%
        dplyr::filter(auc.val.avg == max(auc.val.avg))
      
      output[[i]] <- opt.param
      metrics <- dplyr::bind_rows(output)
      
      # get optimal model per iteration
      opt.model <- ENMeval::eval.models(eval)[[opt.param$tune.args]]
      models[[i]] <- opt.model
      
      # get variable importance for each best model
      var.imp <- ENMeval::eval.variable.importance(eval)[[opt.param$tune.args]]
      contrib[[i]] <- var.imp
      
      # get optimal predictions per iteration
      opt.pred <- ENMeval::eval.predictions(eval)[[opt.param$tune.args]]
      preds[[i]] <- opt.pred
      preds.stack <- raster::stack(preds)
    }
  }
  else if (type == 'type2') {
    for (i in 1:length(bg.list)) {
      
      # make models
      eval <- ENMeval::ENMevaluate(taxon.name = taxon.name, occs = occs, envs = envs, bg = bg.list[[i]], 
                                   tune.args = tune.args, partitions = partitions, partition.settings = partition.settings,
                                   user.grp = user.grp[[i]], doClamp = T, algorithm = 'maxent.jar', parallel = T, parallelType = 'doSNOW')
      
      # get results
      eval.res <- ENMeval::eval.results(eval)
      
      # get optimal parameter combinations
      opt.param <- eval.res %>% dplyr::filter(delta.AICc <= 2) %>%
        dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
        dplyr::filter(auc.val.avg == max(auc.val.avg))
      
      output[[i]] <- opt.param
      metrics <- dplyr::bind_rows(output)
      
      # get optimal model per iteration
      opt.model <- ENMeval::eval.models(eval)[[opt.param$tune.args]]
      models[[i]] <- opt.model
      
      # get variable importance for each best model
      var.imp <- ENMeval::eval.variable.importance(eval)[[opt.param$tune.args]]
      contrib[[i]] <- var.imp
      
      # get optimal predictions per iteration
      opt.pred <- ENMeval::eval.predictions(eval)[[opt.param$tune.args]]
      preds[[i]] <- opt.pred
      preds.stack <- raster::stack(preds)
    }
  }
  return(list(metrics = metrics, models = models, preds = preds.stack, contrib = contrib))
}


### tuning arguments
tune.args <- list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                  rm = seq(0.5,5, by = 0.5))

### bg list
bg.list <- list(bg1_5000, bg1_10000, bg1_15000, bg2_5000, bg2_10000, bg2_15000)


##### ------------------------------------------------------------------------------------------------------------------------------------------------
### O.koreanus == LQ 1.0
# run
o.models_clim <- test_models(taxon.name = 'O.koreanus', occs = o.occs, envs = envs, bg.list = bg.list, tune.args = tune.args,
                             partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(4,4)), type = 'type1') 

# look at metric
print(o.models_clim$metrics)

# look at variable importance
print(o.models_clim$contrib[[2]])

# look at prediction
names(o.models_clim$preds) = c('bg1_5000', 'bg1_10000', 'bg1_15000', 'bg2_5000', 'bg2_10000', 'bg2_15000')
plot(o.models_clim$preds)

# save models
saveRDS(o.models_clim, 'tuning_experiments/output_model_rds/O_koreanus_clim_only_WorldClim.rds')


### K.koreana == LP 5.0
# run
k.models_clim <- test_models(taxon.name = 'K.koreana', occs = k.occs, envs = envs, bg.list = bg.list, tune.args = tune.args,
                             partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(4,4)), type = 'type1')

# look at metric
print(k.models_clim$metrics)

# look at variable importance
print(k.models_clim$contrib[[2]])

# look at prediction
names(k.models_clim$preds) = c('bg1_5000', 'bg1_10000', 'bg1_15000', 'bg2_5000', 'bg2_10000', 'bg2_15000')
plot(k.models_clim$preds)

# save models
saveRDS(k.models_clim, 'tuning_experiments/output_model_rds/K_koreana_clim_only_WorldClim.rds')


#####  Part 8 ::: climate-only model binary ----------------------------------------------------------------------------------------------------------
# calculate thresholds == this function is available in ::: https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/

# ------------------------------------------------------------------------------------------------------------------------
sdm_threshold <- function(sdm, occs, type = "mtp", binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}
# ------------------------------------------------------------------------------------------------------------------------

#### get the list of p10 thresholds
# O. koreanus
o.thresh <- list()

for (i in 1:nlayers(o.models_clim$preds)) {
  thresh <- sdm_threshold(sdm = o.models_clim$preds[[i]], occs = o.occs, type = 'p10', binary = F)
  thresh2 <- raster::minValue(thresh)
  o.thresh[[i]] <- thresh2
  print(o.thresh)
}


# K. koreana
k.thresh <- list()

for (i in 1:nlayers(k.models_clim$preds)) {
  thresh <- sdm_threshold(sdm = k.models_clim$preds[[i]], occs = k.occs, type = 'p10', binary = F)
  thresh2 <- raster::minValue(thresh)
  k.thresh[[i]] <- thresh2
  print(k.thresh)
}


#### automate binary making

# ------------------------------------------------------------------------------------------------------------------------
bin_maker <- function(preds, th) {
  binary.maps <- list()
  
  for (i in 1:nlayers(preds)) {
    bin <- ecospat::ecospat.binary.model(Pred = terra::rast(preds[[i]]), Threshold = th[[i]]) %>% raster::raster()
    binary.maps[[i]] <- bin
    binary.out <- raster::stack(binary.maps)
  }
  return(binary.out)
}
# ------------------------------------------------------------------------------------------------------------------------

#### get binary maps
# O. koreanus
o.bin_clim <- bin_maker(preds = o.models_clim$preds, th = o.thresh)
plot(o.bin_clim)

# K. koreana
k.bin_clim <- bin_maker(preds = k.models_clim$preds, th = k.thresh)
plot(k.bin_clim)


#####  Part 9 ::: plot tuning output ----------------------------------------------------------------------------------------------------------

### O. koreanus continuous
# convert to SpatRaster for layer renaming
o.cont.clim.spat <- terra::rast(o.models_clim$preds)

# recode layer names
names(o.cont.clim.spat) = dplyr::recode(names(o.cont.clim.spat),
                                        'bg1_5000' = 'BG1 (n = 5000)',
                                        'bg1_10000' = 'BG1 (n = 10000)',
                                        'bg1_15000' = 'BG1 (n = 15000)',
                                        'bg2_5000' = 'BG2 (n = 5000)',
                                        'bg2_10000' = 'BG2 (n = 10000)',
                                        'bg2_15000' = 'BG2 (n = 15000)')

# plot
gplot(o.cont.clim.spat) +  
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 1, fill = NA) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/WorldClim tuning/O.koreanus_WorldClim_cont_clim_only.png', width = 30, height = 22, dpi = 800, units = 'cm')


### O. koreanus binary
# convert to SpatRaster for layer renaming
o.bin.clim.spat <- terra::rast(o.bin_clim)

# recode layer names
names(o.bin.clim.spat) = dplyr::recode(names(o.bin.clim.spat),
                                       'bg1_5000' = 'BG1 (n = 5000)',
                                       'bg1_10000' = 'BG1 (n = 10000)',
                                       'bg1_15000' = 'BG1 (n = 15000)',
                                       'bg2_5000' = 'BG2 (n = 5000)',
                                       'bg2_10000' = 'BG2 (n = 10000)',
                                       'bg2_15000' = 'BG2 (n = 15000)')

# plot
gplot(o.bin.clim.spat) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = rev(terrain.colors(1000)),
                       na.value = NA,
                       name = 'Suitability') +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 1, fill = NA) + 
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.position = 'none',
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/WorldClim tuning/O.koreanus_WorldClim_bin_clim_only.png', width = 30, height = 22, dpi = 800, units = 'cm')

#-----------------------------------------------------------------------------------------------------------------------------

### K. koreana continuous
# convert to SpatRaster for layer renaming
k.cont.clim.spat <- terra::rast(k.models_clim$preds)

# recode layer names
names(k.cont.clim.spat) = dplyr::recode(names(k.cont.clim.spat),
                                        'bg1_5000' = 'BG1 (n = 5000)',
                                        'bg1_10000' = 'BG1 (n = 10000)',
                                        'bg1_15000' = 'BG1 (n = 15000)',
                                        'bg2_5000' = 'BG2 (n = 5000)',
                                        'bg2_10000' = 'BG2 (n = 10000)',
                                        'bg2_15000' = 'BG2 (n = 15000)')

# plot
gplot(k.cont.clim.spat) +  
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 1, fill = NA) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/WorldClim tuning/K.koreana_WorldClim_cont_clim_only.png', width = 30, height = 22, dpi = 800, units = 'cm')


### K. koreana binary
# convert to SpatRaster for layer renaming
k.bin.clim.spat <- terra::rast(k.bin_clim)

# recode layer names
names(k.bin.clim.spat) = dplyr::recode(names(k.bin.clim.spat),
                                       'bg1_5000' = 'BG1 (n = 5000)',
                                       'bg1_10000' = 'BG1 (n = 10000)',
                                       'bg1_15000' = 'BG1 (n = 15000)',
                                       'bg2_5000' = 'BG2 (n = 5000)',
                                       'bg2_10000' = 'BG2 (n = 10000)',
                                       'bg2_15000' = 'BG2 (n = 15000)')

# plot
gplot(k.bin.clim.spat) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = rev(terrain.colors(1000)),
                       na.value = NA,
                       name = 'Suitability') +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 1, fill = NA) + 
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.position = 'none',
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/WorldClim tuning/K.koreana_WorldClim_bin_clim_only.png', width = 30, height = 22, dpi = 800, units = 'cm')



#####  Part 10 ::: export everything ----------------------------------------------------------------------------------------------------------
### export model metrics and variable importance
# metrics
print(o.models_clim$metrics)
print(k.models_clim$metrics)

write.csv(o.models_clim$metrics, 'tuning_experiments/metrics/O.koreanus_clim_only_WorldClim_metrics.csv')
write.csv(k.models_clim$metrics, 'tuning_experiments/metrics/K.koreana_clim_only_WorldClim_metrics.csv')

# variable importance
print(o.models_clim$contrib[[2]])
print(k.models_clim$contrib[[2]])

write.csv(o.models_clim$contrib[[2]], 'tuning_experiments/varimp/WorldClim/O.koreanus_clim_only_WorldClim_var.imp.csv')
write.csv(k.models_clim$contrib[[2]], 'tuning_experiments/varimp/WorldClim/K.koreana_clim_only_WorldClim_var.imp.csv')


### export prediction layers 
# O. koreanus cont
for (i in 1:nlayers(o.models_clim$preds)) {
  r <- o.models_clim$preds[[i]]
  file.name <- paste0('tuning_experiments/preds/O.koreanus/WorldClim/cont_clim_only/', names(o.models_clim$preds)[i], '.tif')
  writeRaster(r, file.name, overwrite = T)
}

# O. koreanus bin
for (i in 1:nlayers(o.bin_clim)) {
  r <- o.bin_clim[[i]]
  file.name <- paste0('tuning_experiments/preds/O.koreanus/WorldClim/bin_clim_only/', names(o.bin_clim)[i], '.tif')
  writeRaster(r, file.name, overwrite = T)
}

#-----------------------------------------------------------------------------------------------------------------------------

# K. koreana cont
for (i in 1:nlayers(k.models_clim$preds)) {
  r <- k.models_clim$preds[[i]]
  file.name <- paste0('tuning_experiments/preds/K.koreana/WorldClim/cont_clim_only/', names(k.models_clim$preds)[i], '.tif')
  writeRaster(r, file.name, overwrite = T)
}

# K. koreana bin
for (i in 1:nlayers(k.bin_clim)) {
  r <- k.bin_clim[[i]]
  file.name <- paste0('tuning_experiments/preds/K.koreana/WorldClim/bin_clim_only/', names(k.bin_clim)[i], '.tif')
  writeRaster(r, file.name, overwrite = T)
}

