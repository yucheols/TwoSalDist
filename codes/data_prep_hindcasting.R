#################  hindcasting == use layers acquired from PaleoClim DB
# clear working environment
rm(list = ls(all.names = T))
gc()

# turn off scientific notation
options(scipen = 999)

# load packages
library(raster)
library(ENMeval)
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
# mPWP
mpwp.bio1 <- mpwp[['bio1']]/10
mpwp.bio4 <- mpwp[['bio4']]/10

mpwp <- dropLayer(mpwp, c('bio1', 'bio4'))
mpwp <- raster::stack(mpwp.bio1, mpwp.bio4, mpwp)
print(mpwp)

unitCheck(ref.env = envs, proj.env = mpwp, n = 10000)

# MIS19
# LIG
# LGM
# MH



##### Part 15 ::: fit climate-only model
## load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv') %>% select('long', 'lat')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv') %>% select('long', 'lat')

head(o.occs)
head(k.occs)

## load bg
bg <- read.csv('data/bg/set1/bg1_10000.csv') %>% select('long', 'lat')
head(bg)

## load envs
print(envs)

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

