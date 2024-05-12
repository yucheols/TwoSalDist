##########  plot maps
# clean up working env
rm(list = ls(all.names = T))
gc()

# load packages
library(rasterVis)
library(terra)
library(ggplot2)
library(pals)

### import polygon
poly <- rgdal::readOGR('data/polygons/kor_mer.shp')


#####  Current model

#### O.koreanus
### contunuous
o.clim.wc <- rast('tuning_experiments/preds/O.koreanus/WorldClim/cont_clim_only/bg1_10000.tif')                 # WorldClim clim only
o.full.wc <- rast('tuning_experiments/preds/O.koreanus/WorldClim/cont/O.koreanus_full_model_fixed_parm.tif')    # WorldClim full
o.clim.ch <- rast('tuning_experiments/preds/O.koreanus/CHELSA/cont_clim_only/bg1_10000.tif')                    # CHELSA clim only
o.full.ch <- rast('tuning_experiments/preds/O.koreanus/CHELSA/cont/O.koreanus_full_model_fixed_parm.tif')       # CHELSA full

# combine
o.cont <- c(o.clim.wc, o.full.wc, o.clim.ch, o.full.ch)
names(o.cont) = c('WorldClim', 'WorldClim - full', 'CHELSA', 'CHELSA - full')

# plot
gplot(o.cont) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, ncol = 4, nrow = 1) +
  scale_fill_gradientn(colors =  c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/current/O.koreanus_model_preds.png', width = 20, height = 10, dpi = 800, units = 'cm')


### binary
o.clim.bin.wc <- rast('tuning_experiments/preds/O.koreanus/WorldClim/bin_clim_only/bg1_10000.tif')                 # WorldClim clim only
o.full.bin.wc <- rast('tuning_experiments/preds/O.koreanus/WorldClim/bin/O.koreanus_full_model_fixed_parm.tif')    # WorldClim full
o.clim.bin.ch <- rast('tuning_experiments/preds/O.koreanus/CHELSA/bin_clim_only/bg1_10000.tif')                    # CHELSA clim only
o.full.bin.ch <- rast('tuning_experiments/preds/O.koreanus/CHELSA/bin/O.koreanus_full_model_fixed_parm.tif')       # CHELSA full

# combine
o.bin <- c(o.clim.bin.wc, o.full.bin.wc, o.clim.bin.ch, o.full.bin.ch)
names(o.bin) = c('WorldClim', 'WorldClim - full', 'CHELSA', 'CHELSA - full')

# plot
gplot(o.bin) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, ncol = 4, nrow = 1) +
  scale_fill_gradientn(colors =  rev(terrain.colors(1000)),
                       na.value = NA) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.position = 'none',
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/current/O.koreanus_model_binary.png', width = 20, height = 10, dpi = 800, units = 'cm')

#-------------------------------------------------------------------------------------------------------------------------------------------------------

#### K. koreana
### Continuous
k.clim.wc <- rast('tuning_experiments/preds/K.koreana/WorldClim/cont_clim_only/bg1_10000.tif')                 # WorldClim clim only
k.full.wc <- rast('tuning_experiments/preds/K.koreana/WorldClim/cont/K.koreana_full_model_fixed_parm.tif')    # WorldClim full
k.clim.ch <- rast('tuning_experiments/preds/K.koreana/CHELSA/cont_clim_only/bg1_10000.tif')                    # CHELSA clim only
k.full.ch <- rast('tuning_experiments/preds/K.koreana/CHELSA/cont/K.koreana_full_model_fixed_parm.tif')       # CHELSA full

# combine
k.cont <- c(k.clim.wc, k.full.wc, k.clim.ch, k.full.ch)
names(k.cont) = c('WorldClim', 'WorldClim - full', 'CHELSA', 'CHELSA - full')

# plot
gplot(k.cont) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, ncol = 4, nrow = 1) +
  scale_fill_gradientn(colors =  c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/current/K.koreana_model_preds.png', width = 20, height = 10, dpi = 800, units = 'cm')

#-------------------------------------------------------------------------------------------------------------------------------------------------------

### binary
k.clim.bin.wc <- rast('tuning_experiments/preds/K.koreana/WorldClim/bin_clim_only/bg1_10000.tif')                 # WorldClim clim only
k.full.bin.wc <- rast('tuning_experiments/preds/K.koreana/WorldClim/bin/K.koreana_full_model_fixed_parm.tif')    # WorldClim full
k.clim.bin.ch <- rast('tuning_experiments/preds/K.koreana/CHELSA/bin_clim_only/bg1_10000.tif')                    # CHELSA clim only
k.full.bin.ch <- rast('tuning_experiments/preds/K.koreana/CHELSA/bin/K.koreana_full_model_fixed_parm.tif')       # CHELSA full

# combine
k.bin <- c(k.clim.bin.wc, k.full.bin.wc, k.clim.bin.ch, k.full.bin.ch)
names(k.bin) = c('WorldClim', 'WorldClim - full', 'CHELSA', 'CHELSA - full')

# plot
gplot(k.bin) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, ncol = 4, nrow = 1) +
  scale_fill_gradientn(colors =  rev(terrain.colors(1000)),
                       na.value = NA) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.position = 'none',
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/current/K.koreana_model_binary.png', width = 20, height = 10, dpi = 800, units = 'cm')


#-------------------------------------------------------------------------------------------------------------------------------------------------------

#####  plot hindcasting prediction maps ----------------------------------------------------------------------------------------------------------------

### O.koreanus
# WorldClim based preds
o.hinds.wc <- rast(c(list.files(path = 'hindcast/WorldClim/O_koreanus', pattern = '.tif$', full.names = T)))
sources(o.hinds.wc)
plot(o.hinds.wc[[1]])

names(o.hinds.wc) = c('LGM - WorldClim', 'LIG - WorldClim', 'MH - WorldClim', 'MIS19 - WorldClim', 'mPWP - WorldClim')

# CHELSA based preds
o.hinds.ch <- rast(c(list.files(path = 'hindcast/CHELSA/O_koreanus', pattern = '.tif$', full.names = T)))
o.hinds.ch <- terra::subset(o.hinds.ch, c(2,4,6,8,10))
sources(o.hinds.ch)
plot(o.hinds.ch[[1]])

names(o.hinds.ch) = c('LGM - CHELSA', 'LIG - CHELSA', 'MH - CHELSA', 'MIS19 - CHELSA', 'mPWP - CHELSA')

# combine
o.hinds.comb <- c(o.hinds.wc, o.hinds.ch)
#o.hinds.comb <- terra::subset(o.hinds.comb, c(5,10,4,9,2,7,1,6,3,8))  # for 5 X 2 display
o.hinds.comb <- terra::subset(o.hinds.comb, c(5,4,2,1,3,10,9,7,6,8))  # for 2X 5 display
names(o.hinds.comb)

# plot
gplot(o.hinds.comb) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 2, ncol = 5) +
  scale_fill_gradientn(colors =  c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/hindcast/hindcast_O_koreanus.png', width = 40, height = 15, dpi = 800, units = 'cm')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

### K.koreana
# WorldClim based preds
k.hinds.wc <- rast(c(list.files(path = 'hindcast/WorldClim/K_koreana', pattern = '.tif$', full.names = T)))
sources(k.hinds.wc)
plot(k.hinds.wc[[1]])

names(k.hinds.wc) = c('LGM - WorldClim', 'LIG - WorldClim', 'MH - WorldClim', 'MIS19 - WorldClim', 'mPWP - WorldClim')

# CHELSA based preds
k.hinds.ch <- rast(c(list.files(path = 'hindcast/CHELSA/K_koreana', pattern = '.tif$', full.names = T)))
k.hinds.ch <- terra::subset(k.hinds.ch, c(2,4,6,8,10))
sources(k.hinds.ch)
plot(k.hinds.ch[[1]])

names(k.hinds.ch) = c('LGM - CHELSA', 'LIG - CHELSA', 'MH - CHELSA', 'MIS19 - CHELSA', 'mPWP - CHELSA')

# combine
k.hinds.comb <- c(k.hinds.wc, k.hinds.ch) 
#k.hinds.comb <- terra::subset(k.hinds.comb, c(5,10,4,9,2,7,1,6,3,8))  # for 5 X 2 display
k.hinds.comb <- terra::subset(k.hinds.comb, c(5,4,2,1,3,10,9,7,6,8))  # for 2X 5 display
names(k.hinds.comb)

# plot
gplot(k.hinds.comb) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 2, ncol = 5) +
  scale_fill_gradientn(colors =  c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/hindcast/hindcast_K_koreana.png', width = 40, height = 15, dpi = 800, units = 'cm')


#####  plot MESS ----------------------------------------------------------------------------------------------------------------

##### O.koreanus
## WorldClim
o.mess.wc <- rast(c(list.files(path = 'output_other/WorldClim/clim_only/MESS/O.koreanus', pattern = '.tif$', full.names = T)))
sources(o.mess.wc)
plot(o.mess.wc[[1]])

names(o.mess.wc) = c('LGM - WorldClim', 'LIG - WorldClim', 'MH - WorldClim', 'MIS19 - WorldClim', 'mPWP - WorldClim')
o.mess.wc <- terra::subset(o.mess.wc, c(5,4,2,1,3))

# plot
gplot(o.mess.wc) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 2, ncol = 5) +
  scale_fill_gradientn(colors = rev(as.vector(ocean.thermal(1000))),
                       trans = 'reverse',
                       name = 'MESS',
                       na.value = NA,
                       breaks = c(-1300, -3400),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/hindcast/O.koreanus_MESS_WorldClim.png', width = 40, height = 10, dpi = 800, units = 'cm')


## CHELSA
o.mess.ch <- rast(c(list.files(path = 'output_other/CHELSA/clim_only/MESS/O.koreanus', pattern = '.tif$', full.names = T)))
sources(o.mess.ch)
plot(o.mess.ch[[1]])

names(o.mess.ch) = c('LGM - CHELSA', 'LIG - CHELSA', 'MH - CHELSA', 'MIS19 - CHELSA', 'mPWP - CHELSA')
o.mess.ch <- terra::subset(o.mess.ch, c(5,4,2,1,3))

# plot
gplot(o.mess.ch) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 2, ncol = 5) +
  scale_fill_gradientn(colors = rev(as.vector(ocean.thermal(1000))),
                       trans = 'reverse',
                       name = 'MESS',
                       na.value = NA,
                       breaks = c(-12000, -27500),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/hindcast/O.koreanus_MESS_CHELSA.png', width = 40, height = 10, dpi = 800, units = 'cm')


##### K.koreana
## WorldClim
k.mess.wc <- rast(c(list.files(path = 'output_other/WorldClim/clim_only/MESS/K.koreana', pattern = '.tif$', full.names = T)))
sources(k.mess.wc)
plot(k.mess.wc[[1]])

names(k.mess.wc) = c('LGM - WorldClim', 'LIG - WorldClim', 'MH - WorldClim', 'MIS19 - WorldClim', 'mPWP - WorldClim')
k.mess.wc <- terra::subset(k.mess.wc, c(5,4,2,1,3))

# plot
gplot(k.mess.wc) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 2, ncol = 5) +
  scale_fill_gradientn(colors = rev(as.vector(ocean.thermal(1000))),
                       trans = 'reverse',
                       name = 'MESS',
                       na.value = NA,
                       breaks = c(-2300, -5700),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/hindcast/K.koreana_MESS_WorldClim.png', width = 40, height = 10, dpi = 800, units = 'cm')


## CHELSA
k.mess.ch <- rast(c(list.files(path = 'output_other/CHELSA/clim_only/MESS/K.koreana', pattern = '.tif$', full.names = T)))
sources(k.mess.ch)
plot(k.mess.ch[[1]])

names(k.mess.ch) =  c('LGM - CHELSA', 'LIG - CHELSA', 'MH - CHELSA', 'MIS19 - CHELSA', 'mPWP - CHELSA')
k.mess.ch <- terra::subset(k.mess.ch, c(5,4,2,1,3))

# plot
gplot(k.mess.ch) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 2, ncol = 5) +
  scale_fill_gradientn(colors = rev(as.vector(ocean.thermal(1000))),
                       trans = 'reverse',
                       name = 'MESS',
                       na.value = NA,
                       breaks = c(-23000, -56000),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = poly, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

# save
ggsave('plots/hindcast/K.koreana_MESS_CHELSA.png', width = 40, height = 10, dpi = 800, units = 'cm')
