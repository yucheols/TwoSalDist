# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(raster)
library(dplyr)
library(ENMeval)

##### Part 13 ::: fit climate-only model ---------------------------------------------------------------------------------------------
####  prep data 
## load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv') %>% select('long', 'lat')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv') %>% select('long', 'lat')

head(o.occs)
head(k.occs)

## load bg
# set1
bg1_5000 <- read.csv('data/bg/set1/bg1_5000.csv') %>% select('long', 'lat')
bg1_10000 <- read.csv('data/bg/set1/bg1_10000.csv') %>% select('long', 'lat')
bg1_15000 <- read.csv('data/bg/set1/bg1_15000.csv') %>% select('long', 'lat')

# set2
bg2_5000 <- read.csv('data/bg/set2/bg2_5000.csv') %>% select('long', 'lat')
bg2_10000 <- read.csv('data/bg/set2/bg2_10000.csv') %>% select('long', 'lat')
bg2_15000 <- read.csv('data/bg/set2/bg2_15000.csv') %>% select('long', 'lat')

## load envs
envs <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil', full.names = T))
envs <- raster::stack(subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
print(envs)
plot(envs[[1]])


#### Model tuning
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
### O.koreanus == LQP 3.5
# run
o.models_clim <- test_models(taxon.name = 'O.koreanus', occs = o.occs, envs = envs, bg.list = bg.list, tune.args = tune.args,
                             partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(4,4)), type = 'type1') 

# look at metric
print(o.models_clim$metrics)

# look at variable importance
print(o.models_clim$contrib)

# look at prediction
plot(o.models_clim$preds)

# save models
saveRDS(o.models_clim, 'tuning_experiments/output_model_rds/O_koreanus_clim_only_WorldClim.rds')

# export contribution
write.csv(o.models_clim$contrib, 'tuning_experiments/varimp/WorldCLim/O.koreanus_clim_only_WorldClim_var.imp.csv')

# export metric
write.csv(o.models_clim$metrics, 'tuning_experiments/metrics/O.koreanus_clim_only_WorldClim_metrics.csv')

# export continuous pred
writeRaster(o.models_clim$preds, 'tuning_experiments/preds/O.koreanus/WorldClim/cont/O.koreanus_clim_only.tif', overwrite = T)

##### ------------------------------------------------------------------------------------------------------------------------------------------------


##### ------------------------------------------------------------------------------------------------------------------------------------------------
### K.koreana == Q 0.5
# run
k.models_clim <- test_models(taxon.name = 'K.koreana', occs = k.occs, envs = envs, bg.list = bg.list, tune.args = tune.args,
                             partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(4,4)), type = 'type1')

# look at metric
print(k.models_clim$metrics)

# look at variable importance
print(k.models_clim$contrib)

# look at prediction
plot(k.models_clim$preds)

# save models
saveRDS(k.models_clim, 'tuning_experiments/output_model_rds/K_koreana_clim_only_WorldClim.rds')

# export contribution
write.csv(k.models_clim$contrib, 'tuning_experiments/varimp/WorldClim/K.koreana_clim_only_WorldClim_var.imp.csv')

# export metric
write.csv(k.models_clim$metrics, 'tuning_experiments/metrics/K.koreana_clim_only_WorldClim_metrics.csv')

# export pred
writeRaster(k.models_clim$preds, 'tuning_experiments/preds/K.koreana/WorldClim/cont/K.koreana_clim_only.tif', overwrite = T)

##### ------------------------------------------------------------------------------------------------------------------------------------------------

#####  Part 14 ::: climate-only model binary ----------------------------------------------------------------------------------------------------------
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

##### get thresholds 
# O.koreanus == p10 == 0.3719397 
sdm_threshold(sdm = o.models_clim$preds, occs = o.occs, type = 'p10', binary = F)

# K.koreana == p10 == 0.4518484
sdm_threshold(sdm = k.models_clim$preds, occs = k.occs, type = 'p10', binary = F)


##### get binary
# O.koreanus
o.clim_bin <- ecospat::ecospat.binary.model(Pred = terra::rast(o.models_clim$preds), Threshold = 0.3719397) %>% raster()
plot(o.clim_bin)

writeRaster(o.clim_bin, 'tuning_experiments/preds/O.koreanus/WorldClim/bin/O.koreanus_clim_only_bin.tif', overwrite = T)

# K.koreana 
k.clim_bin <- ecospat::ecospat.binary.model(Pred = terra::rast(k.models_clim$preds), Threshold = 0.4518484) %>% raster()
plot(k.clim_bin)

writeRaster(k.clim_bin, 'tuning_experiments/preds/K.koreana/WorldClim/bin/K.koreana_clim_only_bin.tif', overwrite = T)
