#####  post modeling stuff == model eval, var contributions, resp curves, etc.
#####  this continues from the "clim_only_model_WorldClim.R" workflow
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
o.models_clim <- readRDS('tuning_experiments/output_model_rds/O_koreanus_clim_only_WorldClim.rds')
k.models_clim <- readRDS('tuning_experiments/output_model_rds/K_koreana_clim_only_WorldClim.rds')

glimpse(o.models_clim)
glimpse(k.models_clim)

# check model metrics
print(o.models_clim$metrics)
print(k.models_clim$metrics)

# look at preds
plot(o.models_clim$preds)
plot(k.models_clim$preds)


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
envs <- raster::stack(subset(envs, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))
print(envs)


#####  Part 11 ::: model eval using null models ---------------------------------------------------------------------------------------------
# if user-specified folds were used to fit the models, then you need to provide 'user.eval.type' argument in ENMnulls function

####  O.koreanus null model testing == LQ 1.0
# make ENMeval object as input for ENMnulls
o.e <- ENMevaluate(taxon.name = 'O.koreanus', occs = o.occs[, -1], envs = envs, bg = bg1_10000[, -1], tune.args = list(fc = 'LQ', rm = 1.0), 
                   algorithm = 'maxent.jar', doClamp = T, partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(4,4)))

# test nulls
o.nulls <- ENMnulls(e = o.e, mod.settings =  list(fc = 'LQ', rm = 1.0), eval.stats = c('auc.val', 'auc.diff', 'cbi.val', 'or.10p'), no.iter = 1000)
print(o.nulls@null.results)
print(o.nulls@null.emp.results)

####  K.koreana null model testing == LP 5.0
k.e <- ENMevaluate(taxon.name = 'K.koreana', occs = k.occs[, -1], envs = envs, bg = bg1_10000[, -1], tune.args = list(fc = 'LP', rm = 5.0), 
                   algorithm = 'maxent.jar', doClamp = T, partitions = 'checkerboard2', partition.settings = list(aggregation.factor = c(4,4)))

# test nulls
k.nulls <- ENMnulls(e = k.e, mod.settings =  list(fc = 'LP', rm = 5.0), eval.stats = c('auc.val', 'auc.diff', 'cbi.val', 'or.10p'), no.iter = 1000)
print(k.nulls@null.results)
print(k.nulls@null.emp.results)

####  save null models & results
saveRDS(o.nulls, 'output_nulls/O.koreanus_null_clim_only_WorldClim.rds')
saveRDS(k.nulls, 'output_nulls/K.koreana_null_clim_only_WorldClim.rds')

write.csv(o.nulls@null.results, 'output_nulls/O.koreanus_null_results.csv')
write.csv(o.nulls@null.emp.results, 'output_nulls/O.koreanus_null_summary.csv')

write.csv(k.nulls@null.results, 'output_nulls/K.koreana_null_results.csv')
write.csv(k.nulls@null.emp.results, 'output_nulls/K.koreana_null_summary.csv')


#### plot null model results
# O.koreanus
evalplot.nulls(e.null = o.nulls, stats = c('auc.val', 'cbi.val'), plot.type = 'violin')

# K.koreana
evalplot.nulls(e.null = k.nulls, stats = c('auc.val', 'cbi.val'), plot.type = 'violin')


#####  Part 12 ::: response curves ---------------------------------------------------------------------------------------------
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
o.resp_clim <- respDataPull(sp.name = 'O.koreanus', model = o.models_clim$models[[2]], names.var = names(envs))
print(o.resp_clim)

# K.koreana
k.resp_clim <- respDataPull(sp.name = 'K.koreana', model = k.models_clim$models[[2]], names.var = names(envs))
print(k.resp_clim)


####  plot response curves
# bind
resp_clim <- rbind(o.resp_clim, k.resp_clim)
glimpse(resp_clim)

# reorder vars
resp_clim$var = factor(resp_clim$var, levels = c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15'))
resp_clim$Species = factor(resp_clim$Species, levels = c('O.koreanus', 'K.koreana'))

# recode variable names
resp_clim$var = dplyr::recode_factor(resp_clim$var,
                                     'bio1' = 'Bio1 (°C)', 'bio4' = 'Bio4', 'bio12' = 'Bio12 (mm)', 'bio13' = 'Bio13 (mm)', 
                                     'bio14' = 'Bio14 (mm)', 'bio15' = 'Bio15', 'forest' = 'Forest cover (%)', 'slope' = 'Slope (°)')

# recode species names
resp_clim$Species = dplyr::recode_factor(resp_clim$Species, 'O.koreanus' = 'O. koreanus', 'K.koreana' = 'K. koreana')

# plot
resp_clim %>%
  ggplot(aes(x = x, y = y, group = Species, color = Species)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c('#6495ED', '#ffe600')) +
  facet_wrap(~ var, scales = 'free', nrow = 2, ncol = 3) +
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
ggsave('plots/WorldClim models/response_curves_clim_only.png', width = 30, height = 22, dpi = 800, units = 'cm')
  

#####  Part 13 ::: compare envs values ---------------------------------------------------------------------------------------------
## full envs data
envs_full <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil$', full.names = T))
envs_full <- raster::stack(subset(envs_full, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15', 'forest', 'slope')))
print(envs_full)

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

o.val <- boxdata(sp.name = 'O.koreanus', envs = envs_full, pts = o.occs[, -1])
print(o.val)

k.val <- boxdata(sp.name = 'K.koreana', envs = envs_full, pts = k.occs[, -1])
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
