#####  post modeling stuff == model eval, var contributions, resp curves, etc.
#####  this continues from the "model_WorldClim.R" workflow
getwd()

# clean up working env
rm(list = ls(all.names = T))
gc()

# load packages
library(raster)
library(dplyr)
library(ENMeval)

### load models
o.models <- readRDS('output_model_rds/O_koreanus_model_tuning_WorldClim.rds')
k.models <- readRDS('output_model_rds/K_koreana_model_tuning_WorldClim.rds')

glimpse(o.models)
glimpse(k.models)

### check model metrics
print(o.models$metrics)
print(k.models$metrics)

### look at preds
plot(o.models$preds)
plot(k.models$preds)

#####  Part 10 ::: get variable importance for each sp. ---------------------------------------------------------------------------------------------
# O. koreanus == bg1_10000
print(o.models$contrib[[2]])

# K. koreana == bg1_10000
print(k.models$contrib[[2]])

# export
write.csv(o.models$contrib[[2]], 'data/varimp/O.koreanus_var.imp.csv')
write.csv(k.models$contrib[[2]], 'data/varimp/K.koreana_var_imp.csv')


#####  Part 11 ::: model eval using null models ---------------------------------------------------------------------------------------------
# load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')

# load bg
bg1_10000 <- read.csv('data/bg/set1/bg1_10000.csv')

# load folds
o.folds <- readRDS('data/folds/WorldClim/O.koreanus_folds.rds')
k.folds <- readRDS('data/folds/WorldClim/K.koreana_folds.rds')

# load envs 
envs <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil$', full.names = T))
envs <- raster::stack(subset(envs, c('bio1', 'bio3', 'bio4', 'bio12', 'bio14', 'forest', 'slope')))
print(envs)


####  O.koreanus null model testing == H 2.5
# make ENMeval object as input for ENMnulls
o.e <- ENMevaluate(taxon.name = 'O.koreanus', occs = o.occs[, -1], envs = envs, bg = bg1_10000[, -1], 
                   tune.args = list(fc = 'H', rm = 2.5), algorithm = 'maxent.jar', doClamp = T, partitions = 'user', user.grp = o.folds[[2]])

# test nulls
o.nulls <- ENMnulls(e = o.e, mod.settings =  list(fc = 'H', rm = 2.5), eval.stats = c("auc.val", "auc.diff", "cbi.val", "or.10p"),
                    user.eval.type = 'kspatial', no.iter = 1000)


####  K.koreana null model testing == Q 1.0
k.e <- ENMevaluate(taxon.name = 'K.koreana', occs = k.occs[, -1], envs = envs, bg = bg1_10000[, -1], 
                   tune.args = list(fc = 'Q', rm = 1.0), algorithm = 'maxent.jar', doClamp = T, partitions = 'user', user.grp = k.folds[[2]])

# test nulls
k.nulls <- ENMnulls(e = k.e, mod.settings =  list(fc = 'Q', rm = 1.0), eval.stats = c("auc.val", "auc.diff", "cbi.val", "or.10p"),
                    user.eval.type = 'kspatial', no.iter = 1000)

#####  Part 12 ::: response curves ---------------------------------------------------------------------------------------------
# function to pull out response data

