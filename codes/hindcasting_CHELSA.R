#######  hindcasting from the climate-only ENMs
# clean up working env
rm(list = ls(all.names = T))
gc()

### load packages
library(raster)
library(dplyr)
library(ENMeval)
library(ggplot2)
library(patchwork)

##### Part 16 ::: hindcasting run --------------------------------------------------------------------------------------------------------------------

### load hindcast layers
# mid-Pliocene Warm Period (mPWP)
mpwp <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/mPWP', pattern = '.bil$', full.names = T))
print(mpwp)
plot(mpwp[[1]])

# Marine Isotope Stage 19 (MIS19)
mis <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/MIS19', pattern = '.bil$', full.names = T))
print(mis)
plot(mis[[1]])

# Last Interglacial (LIG)
lig <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/LIG', pattern = '.bil$', full.names = T))
print(lig)
plot(lig[[1]])

# Last Glacial Maximum (LGM)
lgm <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/LGM', pattern = '.bil$', full.names = T))
print(lgm)
plot(lgm[[1]])

# Mid-Holocene (MH)
mh <- raster::stack(list.files(path = 'data/hindcast_layers/processed_CHELSA/MH', pattern = '.bil$', full.names = T))
print(mh)
plot(mh[[1]])


### load climate-only models
# O. koreanus 
#o.models_clim <- readRDS('tuning_experiments/output_model_rds/O_koreanus_clim_only_CHELSA.rds')
#print(o.models_clim)

# K.koreana
#k.models_clim <- readRDS('tuning_experiments/output_model_rds/K_koreana_clim_only_CHELSA.rds')
#print(k.models_clim)

# O. koreanus 
o.models_clim <- readRDS('tuning_experiments/output_model_rds/O_koreanus_clim_only_CHELSA_fixed_bg_params.rds')
print(o.models_clim)

# K.koreana
k.models_clim <- readRDS('tuning_experiments/output_model_rds/K_koreana_clim_only_CHELSA_fixed_bg_params.rds')
print(k.models_clim)


##### model prediction to past climatic conditions
# make function for prediction
model_predictr <- function(model, preds.list, pred.names) {
  require(dismo)
  
  output <- list()
  
  for (i in 1:length(preds.list)) {
    make.pred <- dismo::predict(object = model, x = preds.list[[i]], progress = 'text')
    output[[i]] <- make.pred
  }
  output.pred <- raster::stack(output)
  names(output.pred) = pred.names
  return(output.pred)
}

### O.koreanus
# run
#o.hinds <- model_predictr(model = o.models_clim$models[[1]], 
#                          preds.list = list(mpwp, mis, lig, lgm, mh), 
#                          pred.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'))

# run
o.hinds <- model_predictr(model = o.models_clim@models$fc.LQ_rm.1, 
                          preds.list = list(mpwp, mis, lig, lgm, mh), 
                          pred.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'))

# print & plot
print(o.hinds)
plot(o.hinds)


### K.koreana
# run
#k.hinds <- model_predictr(model = k.models_clim$models[[2]],
#                          preds.list = list(mpwp, mis, lig, lgm, mh),
#                          pred.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'))

# run
k.hinds <- model_predictr(model = k.models_clim@models$fc.LP_rm.5,
                          preds.list = list(mpwp, mis, lig, lgm, mh),
                          pred.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'))

# print & plot
print(k.hinds)
plot(k.hinds)

### export hindcast rasters
# O.koreanus
#for (i in 1:nlayers(o.hinds)) {
#  r <- o.hinds[[i]]
#  name <- paste0('hindcast/CHELSA/O_koreanus/O.koreanus_', names(o.hinds)[i], '.tif')
#  writeRaster(r, filename = name, overwrite = T)
#}

# K.koreana
#for (i in 1:nlayers(k.hinds)) {
#  r <- k.hinds[[i]]
#  name <- paste0('hindcast/CHELSA/K_koreana/K.koreana_', names(k.hinds)[i], '.tif')
#  writeRaster(r, filename = name, overwrite = T)
#}

# O.koreanus
for (i in 1:nlayers(o.hinds)) {
  r <- o.hinds[[i]]
  name <- paste0('hindcast/CHELSA/O_koreanus/O.koreanus_', names(o.hinds)[i], '_fixed_bg_params.tif')
  writeRaster(r, filename = name, overwrite = T)
}

# K.koreana
for (i in 1:nlayers(k.hinds)) {
  r <- k.hinds[[i]]
  name <- paste0('hindcast/CHELSA/K_koreana/K.koreana_', names(k.hinds)[i], '_fixed_bg_params.tif')
  writeRaster(r, filename = name, overwrite = T)
}


##### Part 17 ::: MESS --------------------------------------------------------------------------------------------------------------------
## occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv') %>% dplyr::select(2,3)
k.occs <- read.csv('data/occs/Karsenia_koreana.csv') %>% dplyr::select(2,3)

head(o.occs)
head(k.occs)

## get reference env 
ref.env <- raster::stack(list.files(path = 'data/masked/CHELSA', pattern = '.bil$', full.names = T))
ref.env <- raster::stack(subset(ref.env, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))

## automate MESS
get.mess <- function(proj.env, proj.names, ref.env, occs) {
  require(dplyr)
  require(raster)
  require(dismo)
  
  output <- list()
  ref.env.val <- raster::extract(ref.env, occs) %>% as.data.frame()
  
  for (i in 1:length(proj.names)) {
    mess.calc <- dismo::mess(proj.env[[i]], ref.env.val, full = F)
    output[[i]] <- mess.calc
  }
  stack.mess <- raster::stack(output)
  names(stack.mess) = proj.names
  
  return(stack.mess)
}

## get MESS 
# O.koreanus
o.mess <- get.mess(proj.env = list(mpwp, mis, lig, lgm, mh), proj.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'), ref.env = ref.env, occs = o.occs)
plot(o.mess)

# K.koreana
k.mess <- get.mess(proj.env = list(mpwp, mis, lig, lgm, mh), proj.names = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH'), ref.env = ref.env, occs = k.occs)
plot(k.mess)

# export MESS
writeRaster(o.mess, 'output_other/CHELSA/clim_only/O.koreanus_clim_only_MESS.tif', overwrite = T)
writeRaster(k.mess, 'output_other/CHELSA/clim_only/K.koreana_clim_only_MESS.tif', overwrite = T)


##### Part 18 :::  compare env values between time ----------------------------------------------------------------------
## use bio1 (most important var for O.koreanus) and bio13 (most important var for K.koreana)

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


##### O.koreanus

### get data
# mPWP
o.dat.mpwp <- boxdata(sp.name = 'O.koreanus', envs = mpwp[[c('bio1', 'bio13')]], pts = o.occs)
o.dat.mpwp$time = 'mPWP'
o.dat.mpwp[1:187, 1] <- o.dat.mpwp[1:187, 1]/10
head(o.dat.mpwp)

# MIS19
o.dat.mis <- boxdata(sp.name = 'O.koreanus', envs = mis[[c('bio1', 'bio13')]], pts = o.occs)
o.dat.mis$time = 'MIS19'
o.dat.mis[1:187, 1] <- o.dat.mis[1:187, 1]/10
head(o.dat.mis)

# LIG
o.dat.lig <- boxdata(sp.name = 'O.koreanus', envs = lig[[c('bio1', 'bio13')]], pts = o.occs)
o.dat.lig$time = 'LIG'
o.dat.lig[1:187, 1] <- o.dat.lig[1:187, 1]/10
head(o.dat.lig)

# LGM
o.dat.lgm <- boxdata(sp.name = 'O.koreanus', envs = lgm[[c('bio1', 'bio13')]], pts = o.occs)
o.dat.lgm$time = 'LGM'
o.dat.lgm[1:187, 1] <- o.dat.lgm[1:187, 1]/10
head(o.dat.lgm)

# MH
o.dat.mh <- boxdata(sp.name = 'O.koreanus', envs = mh[[c('bio1', 'bio13')]], pts = o.occs)
o.dat.mh$time = 'MH'
o.dat.mh[1:187, 1] <- o.dat.mh[1:187, 1]/10
head(o.dat.mh)

# current
o.dat.cur <- boxdata(sp.name = 'O.koreanus', envs = ref.env[[c('bio1', 'bio13')]], pts = o.occs)
o.dat.cur$time = 'Current'
o.dat.cur[1:187, 1] <- o.dat.cur[1:187, 1]/10
head(o.dat.cur)

## format data for plotting
o.dat.bind <- rbind(o.dat.mpwp, o.dat.mis, o.dat.lig, o.dat.lgm, o.dat.mh, o.dat.cur)
head(o.dat.bind)

o.dat.bind$time <- factor(o.dat.bind$time, levels = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'))
o.dat.bind$var <- factor(o.dat.bind$var, levels = c('bio1', 'bio13')) 
o.dat.bind$var <- dplyr::recode(o.dat.bind$var, 'bio1' = 'Bio1 (°C)', 'bio13' = 'Bio13 (mm)')

## plot
o.env.plot <- o.dat.bind %>%
  ggplot(aes(x = time, y = val)) +
  geom_boxplot(fill = '#6495ED', color = '#6495ED', linewidth = 1.0, alpha = 0.4, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(aes(color = var), position = position_jitterdodge(0.4), alpha = 0.4) +
  scale_color_manual(values = rep('#6495ED', 2)) +
  facet_wrap(~ var, scales = 'free') +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = 'none')

##### K.koreana

### get data
# mPWP
#k.dat.mpwp <- boxdata(sp.name = 'K.koreana', envs = mpwp[[c('bio4', 'bio12')]], pts = k.occs)
#k.dat.mpwp$time = 'mPWP'
#head(k.dat.mpwp)

# MIS19
#k.dat.mis <- boxdata(sp.name = 'K.koreana', envs = mis[[c('bio4', 'bio12')]], pts = k.occs)
#k.dat.mis$time = 'MIS19'
#head(k.dat.mis)

# LIG
#k.dat.lig <- boxdata(sp.name = 'K.koreana', envs = lig[[c('bio4', 'bio12')]], pts = k.occs)
#k.dat.lig$time = 'LIG'
#head(k.dat.lig)

# LGM
#k.dat.lgm <- boxdata(sp.name = 'K.koreana', envs = lgm[[c('bio4', 'bio12')]], pts = k.occs)
#k.dat.lgm$time = 'LGM'
#head(k.dat.lgm)

# MH
#k.dat.mh <- boxdata(sp.name = 'K.koreana', envs = mh[[c('bio4', 'bio12')]], pts = k.occs)
#k.dat.mh$time = 'MH'
#head(k.dat.mh)

# current
#k.dat.cur <- boxdata(sp.name = 'K.koreana', envs = ref.env[[c('bio4', 'bio12')]], pts = k.occs)
#k.dat.cur$time = 'Current'
#head(k.dat.cur)

## format data for plotting
#k.dat.bind <- rbind(k.dat.mpwp, k.dat.mis, k.dat.lig, k.dat.lgm, k.dat.mh, k.dat.cur)
#head(k.dat.bind)

#k.dat.bind$time <- factor(k.dat.bind$time, levels = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'))
#k.dat.bind$var <- factor(k.dat.bind$var, levels = c('bio4', 'bio12'))
#k.dat.bind$var <- dplyr::recode(k.dat.bind$var, 'bio4' = 'Bio4', 'bio12' = 'Bio12 (mm)')

### get data
# mPWP
k.dat.mpwp <- boxdata(sp.name = 'K.koreana', envs = mpwp[[c('bio4', 'bio14')]], pts = k.occs)
k.dat.mpwp$time = 'mPWP'
k.dat.mpwp[1:137, 1] <- k.dat.mpwp[1:137, 1]/10
head(k.dat.mpwp)

# MIS19
k.dat.mis <- boxdata(sp.name = 'K.koreana', envs = mis[[c('bio4', 'bio14')]], pts = k.occs)
k.dat.mis$time = 'MIS19'
k.dat.mis[1:137, 1] <- k.dat.mis[1:137, 1]/10
head(k.dat.mis)

# LIG
k.dat.lig <- boxdata(sp.name = 'K.koreana', envs = lig[[c('bio4', 'bio14')]], pts = k.occs)
k.dat.lig$time = 'LIG'
k.dat.lig[1:137, 1] <- k.dat.lig[1:137, 1]/10
head(k.dat.lig)

# LGM
k.dat.lgm <- boxdata(sp.name = 'K.koreana', envs = lgm[[c('bio4', 'bio14')]], pts = k.occs)
k.dat.lgm$time = 'LGM'
k.dat.lgm[1:137, 1] <- k.dat.lgm[1:137, 1]/10
head(k.dat.lgm)

# MH
k.dat.mh <- boxdata(sp.name = 'K.koreana', envs = mh[[c('bio4', 'bio14')]], pts = k.occs)
k.dat.mh$time = 'MH'
k.dat.mh[1:137, 1] <- k.dat.mh[1:137, 1]/10
head(k.dat.mh)

# current
k.dat.cur <- boxdata(sp.name = 'K.koreana', envs = ref.env[[c('bio4', 'bio14')]], pts = k.occs)
k.dat.cur$time = 'Current'
k.dat.cur[1:137, 1] <- k.dat.cur[1:137, 1]/10
head(k.dat.cur)

## format data for plotting
k.dat.bind <- rbind(k.dat.mpwp, k.dat.mis, k.dat.lig, k.dat.lgm, k.dat.mh, k.dat.cur)
head(k.dat.bind)

k.dat.bind$time <- factor(k.dat.bind$time, levels = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'))
k.dat.bind$var <- factor(k.dat.bind$var, levels = c('bio4', 'bio14'))
k.dat.bind$var <- dplyr::recode(k.dat.bind$var, 'bio4' = 'Bio4', 'bio14' = 'Bio14 (mm)')

## plot
k.env.plot <- k.dat.bind %>%
  ggplot(aes(x = time, y = val)) +
  geom_boxplot(fill = '#ffe600', color = '#ffe600', linewidth = 1.0, alpha = 0.4, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(aes(color = var), position = position_jitterdodge(0.4), alpha = 0.4) +
  scale_color_manual(values = rep('#ffe600', 2)) +
  facet_wrap(~ var, scales = 'free') +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = 'none')


## combine plots for the two species into one == using the patchwork library
o.env.plot + k.env.plot +
  plot_layout(nrow = 2)

## export   
#ggsave('plots/CHELSA models/env_values_thru_time.png', width = 20, height = 25, dpi = 800, units = 'cm')
ggsave('plots/CHELSA models/env_values_thru_time_fixed_bg_params.png', width = 20, height = 25, dpi = 800, units = 'cm')


##### Part 19 :::  optional == visualize climate trends through time // from mPWP to current --------------------------------------------------------------
##### annual mean temp and annual precip

## function to automate data formatting
env.val.get <- function(env.list, var, time, type) {
  require(dplyr)
  
  output <- list()
  
  if (type == 'max'){
    for (i in 1:length(env.list)) {
      get.df <- env.list[[i]][[var]] %>% as.data.frame() %>% na.omit()
      get.df.max <- max(get.df[[var]])
      df.max <- data.frame(val = get.df.max, time = time[[i]])
      
      output[[i]] <- df.max
      output[[i]]$type = type
      output[[i]]$var = var
    }
  }
  else if (type == 'mean') {
    for (i in 1:length(env.list)) {
      get.df <- env.list[[i]][[var]] %>% as.data.frame() %>% na.omit()
      get.df.mean <- mean(get.df[[var]])
      df.mean <- data.frame(val = get.df.mean, time = time[[i]])
      
      output[[i]] <- df.mean
      output[[i]]$type = type
      output[[i]]$var = var
    }
  }
  else if (type == 'min') {
    for (i in 1:length(env.list)) {
      get.df <- env.list[[i]][[var]] %>% as.data.frame() %>% na.omit()
      get.df.min <- min(get.df[[var]])
      df.min <- data.frame(val = get.df.min, time = time[[i]])
      
      output[[i]] <- df.min
      output[[i]]$type = type
      output[[i]]$var = var
    }
  }
  output.all <- dplyr::bind_rows(output)
  return(output.all)
}


### load current values
print(ref.env)
plot(ref.env[[1]])

### bio1 ::: Annual Mean Temp
# min
min.vals.tmp <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, ref.env), var = 'bio1', 
                            time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'min')

min.vals.tmp$val <- min.vals.tmp$val/10
print(min.vals.tmp)

# mean
mean.vals.tmp <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, ref.env), var = 'bio1', 
                             time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'mean')

mean.vals.tmp$val <- mean.vals.tmp$val/10
print(mean.vals.tmp)

# max
max.vals.tmp <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, ref.env), var = 'bio1', 
                            time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'max')

max.vals.tmp$val <- max.vals.tmp$val/10
print(max.vals.tmp)


### bio12 ::: Annual Precipitation
# min
min.vals.pr <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, ref.env), var = 'bio12', 
                           time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'min')

print(min.vals.pr)

# mean
mean.vals.pr <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, ref.env), var = 'bio12', 
                            time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'mean')

print(mean.vals.pr)

# max
max.vals.pr <- env.val.get(env.list = list(mpwp, mis, lig, lgm, mh, ref.env), var = 'bio12', 
                           time = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'), type = 'max')

print(max.vals.pr)


### bind
env.vals <- rbind(min.vals.tmp, mean.vals.tmp, max.vals.tmp, min.vals.pr, mean.vals.pr, max.vals.pr)
print(env.vals)

### reorder variable levels 
env.vals$time = factor(env.vals$time, levels = c('mPWP', 'MIS19', 'LIG', 'LGM', 'MH', 'Current'))
env.vals$type = factor(env.vals$type, levels = c('min', 'mean', 'max'))

### rename variables
env.vals$var <- dplyr::recode_factor(env.vals$var, 'bio1' = 'Bio1 (°C)', 'bio12' = 'Bio12 (mm)')
env.vals$type <- dplyr::recode_factor(env.vals$type, 'min' = 'Min', 'mean' = 'Mean', 'max' = 'Max')

### plot
env.vals %>%
  ggplot(aes(x = time, y = val, group = type, color = type)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ var, scales = 'free') +
  xlab('Time period') + ylab('Values') +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = 'top')

### save
ggsave('plots/CHELSA models/clim_trends_thru_time.png', width = 20, height = 25, dpi = 800, units = 'cm')
