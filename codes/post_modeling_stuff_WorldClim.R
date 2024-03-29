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
o.models <- readRDS('tuning_experiments/output_model_rds/O_koreanus_model_tuning_WorldClim.rds')
k.models <- readRDS('tuning_experiments/output_model_rds/K_koreana_model_tuning_WorldClim.rds')

glimpse(o.models)
glimpse(k.models)

# check model metrics
print(o.models$metrics)
print(k.models$metrics)

# look at preds
plot(o.models$preds)
plot(k.models$preds)


### load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')

### load bg
bg1_10000 <- read.csv('data/bg/set1/bg1_10000.csv')

### load folds == only needed if user specified folds were used to make the models
#o.folds <- readRDS('data/folds/WorldClim/O.koreanus_user_folds.rds')
#k.folds <- readRDS('data/folds/WorldClim/K.koreana_folds.rds')

### load envs 
envs <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil$', full.names = T))
envs <- raster::stack(subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15', 'forest', 'slope')))
print(envs)


#####  Part 9 ::: get variable importance for each sp. ---------------------------------------------------------------------------------------------
# O. koreanus == bg1_10000
print(o.models$contrib[[2]])

# K. koreana == bg1_10000
print(k.models$contrib[[2]])

# export
write.csv(o.models$contrib[[2]], 'tuning_experiments/varimp/WorldClim/O.koreanus_WorldClim_var.imp.csv')
write.csv(k.models$contrib[[2]], 'tuning_experiments/varimp/WorldClim/K.koreana_WorldClim_var_imp.csv')


#####  Part 10 ::: model eval using null models ---------------------------------------------------------------------------------------------
# if user-specified folds were used to fit the models, then you need to provide 'user.eval.type' argument in ENMnulls function

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
saveRDS(o.nulls, 'output_nulls/O.koreanus_null_WorldClim.rds')
saveRDS(k.nulls, 'output_nulls/K.koreana_null_WorldClim.rds')


#### plot null model results
# O.koreanus
evalplot.nulls(e.null = o.nulls, stats = c('auc.val', 'cbi.val'), plot.type = 'violin')

# K.koreana
evalplot.nulls(e.null = k.nulls, stats = c('auc.val', 'cbi.val'), plot.type = 'violin')


#####  Part 11 ::: response curves ---------------------------------------------------------------------------------------------
# function to pull out response data
respDataPull <- function(sp.name, model, names.var) {
  require(dplyr)
  
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

# recode variable names
resp$var = dplyr::recode_factor(resp$var,
                                'bio1' = 'Bio1 (°C)', 'bio4' = 'Bio4', 'bio12' = 'Bio12 (mm)', 'bio13' = 'Bio13 (mm)', 
                                'bio14' = 'Bio14 (mm)', 'bio15' = 'Bio15', 'forest' = 'Forest cover (%)', 'slope' = 'Slope (°)')

# recode species names
resp$Species = dplyr::recode_factor(resp$Species, 'O.koreanus' = 'O. koreanus', 'K.koreana' = 'K. koreana')

# plot
resp %>%
  ggplot(aes(x = x, y = y, group = Species, color = Species)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c('#6495ED', '#ffe600')) +
  facet_wrap(~ var, scales = 'free', nrow = 2, ncol = 4) +
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

# save plot
ggsave('plots/WorldClim models/response_curves.png', width = 30, height = 22, dpi = 800, units = 'cm')
  

#####  Part 12 ::: compare envs values ---------------------------------------------------------------------------------------------
## extract envs values
print(envs)

## O. koreanus
# function to format data for box plot
boxdata <- function(sp.name, envs, pts) {
  require(dplyr)
  
  output <- list()
  
  for (i in 1:length(names(envs))) {
    val <- raster::extract(envs[[i]], pts) %>% as.data.frame()
    val$var = names(envs)[i]
    val$species = sp.name
    colnames(val) = c('val', 'var', 'Species')
    output[[i]] <- val
  }
  output <- dplyr::bind_rows(output)
  return(output)
}

o.val <- boxdata(sp.name = 'O.koreanus', envs = envs, pts = o.occs[, -1])
print(o.val)

k.val <- boxdata(sp.name = 'K.koreana', envs = envs, pts = k.occs[, -1])
print(k.val)

vals <- rbind(o.val, k.val)
head(vals)

## reorder species plotting order
vals$Species = factor(vals$Species, levels = c('O.koreanus', 'K.koreana'))

## recode variable names
vals$var = dplyr::recode_factor(vals$var, 
                                'bio1' = 'Bio1 (°C)', 'bio4' = 'Bio4', 'bio12' = 'Bio12 (mm)', 'bio13' = 'Bio13 (mm)', 
                                'bio14' = 'Bio14 (mm)', 'bio15' = 'Bio15', 'forest' = 'Forest cover (%)', 'slope' = 'Slope (°)')

## recode species names
vals$Species = dplyr::recode_factor(vals$Species, 'O.koreanus' = 'O. koreanus', 'K.koreana' = 'K. koreana')

## plot box
vals %>%
  ggplot(aes(x = var, y = val, fill = Species, color = Species)) +
  geom_boxplot(linewidth = 1.0, alpha = 0.4, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(0.4), alpha = 0.4) +
  facet_wrap(~ var, scale = 'free', nrow = 2, ncol = 4) +
  scale_fill_manual(values = c('#6495ED', '#ffe600')) +
  scale_color_manual(values = c('#6495ED', '#ffe600')) +
  xlab('Variable') + ylab('Value') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14, face = 'italic'),
        legend.position = 'top')

## save
ggsave('plots/WorldClim models/env_values_boxplot.png', width = 20, height = 25, dpi = 800, units = 'cm')
