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
library(ggplot2)

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
# if user-specified folds were used to fit the models, then you need to provide 'user.eval.type' argument in ENMnulls function

# load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')

# load bg
bg1_10000 <- read.csv('data/bg/set1/bg1_10000.csv')

# load folds == only needed if user specified folds were used to make the models
#o.folds <- readRDS('data/folds/WorldClim/O.koreanus_user_folds.rds')
#k.folds <- readRDS('data/folds/WorldClim/K.koreana_folds.rds')

# load envs 
envs <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil$', full.names = T))
envs <- raster::stack(subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15', 'forest', 'slope')))
print(envs)


####  O.koreanus null model testing == LP 1.5
# make ENMeval object as input for ENMnulls
o.e <- ENMevaluate(taxon.name = 'O.koreanus', occs = o.occs[, -1], envs = envs, bg = bg1_10000[, -1], tune.args = list(fc = 'LP', rm = 1.5), 
                   algorithm = 'maxent.jar', doClamp = T, partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(7,7)))

# test nulls
o.nulls <- ENMnulls(e = o.e, mod.settings =  list(fc = 'LP', rm = 1.5), eval.stats = c('auc.val', 'auc.diff', 'cbi.val', 'or.10p'), no.iter = 1000)


####  K.koreana null model testing == Q 4.5
k.e <- ENMevaluate(taxon.name = 'K.koreana', occs = k.occs[, -1], envs = envs, bg = bg1_10000[, -1], tune.args = list(fc = 'Q', rm = 4.5), 
                   algorithm = 'maxent.jar', doClamp = T, partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(7,7)))

# test nulls
k.nulls <- ENMnulls(e = k.e, mod.settings =  list(fc = 'Q', rm = 4.5), eval.stats = c('auc.val', 'auc.diff', 'cbi.val', 'or.10p'), no.iter = 1000)


####  save null models
#saveRDS(o.nulls, 'output_nulls/O.koreanus_null.rds')
#saveRDS(k.nulls, 'output_nulls/K.koreana_null.rds')


#### plot null model results
# O.koreanus
evalplot.nulls(e.null = o.nulls, stats = c('auc.val', 'cbi.val'), plot.type = 'violin')

# K.koreana
evalplot.nulls(e.null = k.nulls, stats = c('auc.val', 'cbi.val'), plot.type = 'violin')


#####  Part 12 ::: response curves ---------------------------------------------------------------------------------------------
# function to pull out response data
respDataPull <- function(sp.name, model, names.var) {
  resp.data <- list()
  
  for (i in 1:length(names.var)) {
    resp <- as.data.frame(dismo::response(x = model, var = names.var[[i]]))
    colnames(resp) = c('x','y')
    resp$var = names.var[[i]]
    resp.data[[i]] <- resp
    resp.data.out <- dplyr::bind_rows(resp.data)
    resp.data.out$Species = sp.name
  }
  return(resp.data.out)
}

#####  pull out the data
# O.koreanus
o.resp <- respDataPull(sp.name = 'O.koreanus', model = o.models$models[[2]], names.var = names(envs))
print(o.resp)

# K.koreana
k.resp <- respDataPull(sp.name = 'K.koreana', model = k.models$models[[2]], names.var = names(envs))
print(k.resp)


####  plot response curves
# bind
resp <- rbind(o.resp, k.resp)
glimpse(resp)

# reorder vars
resp$var = factor(resp$var, levels = c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15', 'forest', 'slope'))
resp$Species = factor(resp$Species, levels = c('O.koreanus', 'K.koreana'))

# plot
resp %>%
  ggplot(aes(x = x, y = y, group = Species, color = Species)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c('#6495ED', '#FDEF3B')) +
  facet_wrap(~ var, scales = 'free') +
  xlab('Value') + ylab('Suitability') +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14, face = 'italic'),
        legend.position = 'top')
  

#####  Part 13 ::: compare envs values ---------------------------------------------------------------------------------------------
## extract envs values
print(envs)

o.val <- raster::extract(envs, o.occs[, -1]) %>% as.data.frame()
o.val$Species = 'O.koreanus'

k.val <- raster::extract(envs, k.occs[, -1]) %>% as.data.frame()
k.val$Species = 'K.koreana'

vals <- rbind(o.val, k.val)
print(vals)

## plot box
