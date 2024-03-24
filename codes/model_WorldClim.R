##### Data prep for current ENM with WorldClim data 
# set seed 
set.seed(111)

# load packages
library(ENMeval)
library(dplyr)
library(raster)
library(rasterVis)
library(ggplot2)

# check data
print(envs)

head(o.occs)
head(k.occs)

head(bg1_5000)
head(bg2_5000)

#glimpse(o.folds)
#glimpse(k.folds)

#####  Part 6 ::: model testing  ---------------------------------------------------------------------------------------------
# automate model tuning 
# type 1 == minimum or.10p.avg as primary criterion // type 2 == delta.AICc <= 2 as primary criterion 
test_models <- function(taxon.name, occs, envs, bg.list, tune.args, partitions, partition.settings = NULL, user.grp = NULL, type) {
  require(dplyr)
  require(ENMeval)
  require(raster)
  
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


### prep inputs
bg.list <- list(bg1_5000[, c('long', 'lat')], bg1_10000[, c('long', 'lat')], bg1_15000[, c('long', 'lat')],
                bg2_5000[, c('long', 'lat')], bg2_10000[, c('long', 'lat')], bg2_15000[, c('long', 'lat')])

tune.args <- list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                  rm = seq(0.5,5, by = 0.5))


### O. koreanus model testing run
# run
o.models <- test_models(taxon.name = 'O.koreanus', occs = o.occs[, c(2,3)], envs = envs, bg.list = bg.list, tune.args = tune.args, 
                        partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(7,7)), type = 'type1')

# look at results
print(o.models$metrics)

# look at predictions
names(o.models$preds) = c('bg1_5000', 'bg1_10000', 'bg1_15000', 'bg2_5000', 'bg2_10000', 'bg2_15000')
plot(o.models$preds)

# save output as .rds for later use
saveRDS(o.models, 'output_model_rds/O_koreanus_model_tuning_WorldClim.rds')


### K. koreana model testing run
k.models <- test_models(taxon.name = 'K.koreana', occs = k.occs[, c(2,3)], envs = envs, bg.list = bg.list, tune.args = tune.args, 
                        partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(7,7)), type = 'type1')

# look at results
print(k.models$metrics)

# look at predictions
names(k.models$preds) = c('bg1_5000', 'bg1_10000', 'bg1_15000', 'bg2_5000', 'bg2_10000', 'bg2_15000')
plot(k.models$preds)

# save output as .rds for later use
saveRDS(k.models, 'output_model_rds/K_koreana_model_tuning_WorldClim.rds')


#####  Part 7 ::: look at binary  ---------------------------------------------------------------------------------------------
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

for (i in 1:nlayers(o.models$preds)) {
  thresh <- sdm_threshold(sdm = o.models$preds[[i]], occs = o.occs[, c(2,3)], type = 'p10', binary = F)
  thresh2 <- raster::minValue(thresh)
  o.thresh[[i]] <- thresh2
  print(o.thresh)
}


# K. koreana
k.thresh <- list()

for (i in 1:nlayers(k.models$preds)) {
  thresh <- sdm_threshold(sdm = k.models$preds[[i]], occs = k.occs[, c(2,3)], type = 'p10', binary = F)
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
o.bin <- bin_maker(preds = o.models$preds, th = o.thresh)
plot(o.bin)

# K. koreana
k.bin <- bin_maker(preds = k.models$preds, th = k.thresh)
plot(k.bin)


#####  Part 8 ::: plot tuning outputs  ---------------------------------------------------------------------------------------------

### O. koreanus continuous
# plot
gplot(o.models$preds) +  
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
ggsave('plots/WorldClim tuning/O.koreanus_WorldClim_cont.png', width = 30, height = 22, dpi = 800, units = 'cm')


### O. koreanus binary
# plot
gplot(o.bin) + 
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
ggsave('plots/WorldClim tuning/O.koreanus_WorldClim_bin.png', width = 30, height = 22, dpi = 800, units = 'cm')

#-----------------------------------------------------------------------------------------------------------------------------

### K. koreana continuous
# plot
gplot(k.models$preds) +  
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
ggsave('plots/WorldClim tuning/K.koreana_WorldClim_cont.png', width = 30, height = 22, dpi = 800, units = 'cm')


### K. koreana binary
# plot
gplot(k.bin) + 
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
ggsave('plots/WorldClim tuning/K.koreana_WorldClim_bin.png', width = 30, height = 22, dpi = 800, units = 'cm')

#-----------------------------------------------------------------------------------------------------------------------------

##### export model metrics and prediction layers
### metrics
print(o.models$metrics)
print(k.models$metrics)

write.csv(o.models$metrics, 'tuning_experiments/metrics/O.koreanus_WorldClim_metrics.csv')
write.csv(k.models$metrics, 'tuning_experiments/metrics/K.koreana_WorldClim_metrics.csv')

### prediction layers
# O. koreanus cont
for (i in 1:nlayers(o.models$preds)) {
  r <- o.models$preds[[i]]
  file.name <- paste0('tuning_experiments/preds/O.koreanus/WorldClim/cont/', names(o.models$preds)[i], '.tif')
  writeRaster(r, file.name, overwrite = T)
}

# O. koreanus bin
for (i in 1:nlayers(o.bin)) {
  r <- o.bin[[i]]
  file.name <- paste0('tuning_experiments/preds/O.koreanus/WorldClim/bin/', names(o.bin)[i], '.tif')
  writeRaster(r, file.name, overwrite = T)
}

#-----------------------------------------------------------------------------------------------------------------------------

# K. koreana cont
for (i in 1:nlayers(k.models$preds)) {
  r <- k.models$preds[[i]]
  file.name <- paste0('tuning_experiments/preds/K.koreana/WorldClim/cont/', names(k.models$preds)[i], '.tif')
  writeRaster(r, file.name, overwrite = T)
}

# K. koreana bin
for (i in 1:nlayers(k.bin)) {
  r <- k.bin[[i]]
  file.name <- paste0('tuning_experiments/preds/K.koreana/WorldClim/bin/', names(k.bin)[i], '.tif')
  writeRaster(r, file.name, overwrite = T)
}

