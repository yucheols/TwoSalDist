#####  sptially visulize WorldClim vs. CHELSA prediction differences
#####  use absolute difference to visualize where the models disagree regardless of direction

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(terra)
library(rgdal) 
library(rasterVis)
library(ggplot2)
library(pals)

# load korea polygon
kor <- readOGR('data/polygons/kor_mer.shp')


#####  part 1 ::: O. koreanus  ----------

###  current
# load predictions
o.wc.full <- rast('tuning_experiments/preds/O.koreanus/WorldClim/cont/bg1_10000.tif')
plot(o.wc.full)

o.ch.full <- rast('tuning_experiments/preds/O.koreanus/CHELSA/cont/bg1_10000.tif') 
plot(o.ch.full)

o.wc.clim <- rast('tuning_experiments/preds/O.koreanus/WorldClim/cont_clim_only/bg1_10000.tif')
plot(o.wc.clim)

o.ch.clim <- rast('tuning_experiments/preds/O.koreanus/CHELSA/cont_clim_only/bg1_10000.tif')
plot(o.ch.clim)

# visualize deviations
o.full.dev <- abs(o.wc.full - o.ch.full)
plot(o.full.dev)

o.clim.dev <- abs(o.wc.clim - o.ch.clim)
plot(o.clim.dev)


###  hindcast
# load WorldClim-based predictions
o.wc.hind <- rast(list.files(path = 'hindcast/WorldClim/O_koreanus/', pattern = '.tif$', full.names = T))
sources(o.wc.hind)

names(o.wc.hind) = c('LGM', 'LIG', 'MH', 'MIS19', 'mPWP')
o.wc.hind <- o.wc.hind[[c(5,4,2,1,3)]]
plot(o.wc.hind)

# load CHELSA-based predictions 
o.ch.hind <- rast(list.files(path = 'hindcast/CHELSA/O_koreanus/hindcast_fixed_bg_params/', pattern = '.tif$', full.names = T))
sources(o.ch.hind)

names(o.ch.hind) = c('LGM', 'LIG', 'MH', 'MIS19', 'mPWP')
o.ch.hind <- o.ch.hind[[c(5,4,2,1,3)]]
plot(o.ch.hind)

# visualize deviations
o.hind.dev <- abs(o.wc.hind - o.ch.hind)
print(o.hind.dev)
plot(o.hind.dev)


#####  part 2 ::: K. koreana  ----------

###  current
# load predictions
k.wc.full <- rast('tuning_experiments/preds/K.koreana/WorldClim/cont/bg1_10000.tif')
plot(k.wc.full)

k.ch.full <- rast('tuning_experiments/preds/K.koreana/CHELSA/cont/bg1_10000.tif')
plot(k.ch.full)

k.wc.clim <- rast('tuning_experiments/preds/K.koreana/WorldClim/cont_clim_only/bg1_10000.tif')
plot(k.wc.clim)

k.ch.clim <- rast('tuning_experiments/preds/K.koreana/CHELSA/cont/K.koreana_full_model_fixed_parm.tif')
plot(k.ch.clim)

# visualize deviations
k.full.dev <- abs(k.wc.full - k.ch.full)
plot(k.full.dev)

k.clim.dev <- abs(k.wc.clim - k.ch.clim)
plot(k.clim.dev)


###  hindcast
# load WorldClim-based predictions
k.wc.hind <- rast(list.files(path = 'hindcast/WorldClim/K_koreana/', pattern = '.tif$', full.names = T))
sources(k.wc.hind)

names(k.wc.hind) = c('LGM', 'LIG', 'MH', 'MIS19', 'mPWP')
k.wc.hind <- k.wc.hind[[c(5,4,2,1,3)]]
plot(k.wc.hind)

# load CHELSA-based predictions 
k.ch.hind <- rast(list.files('hindcast/CHELSA/K_koreana/hindcast_fixed_bg_params/', pattern = '.tif$', full.names = T))
sources(k.ch.hind)

names(k.ch.hind) = c('LGM', 'LIG', 'MH', 'MIS19', 'mPWP')
k.ch.hind <- k.ch.hind[[c(5,4,2,1,3)]]
plot(k.ch.hind)

# visualize deviations
k.hind.dev <- abs(k.wc.hind - k.ch.hind)
plot(k.hind.dev)


#####  part 3 ::: plot  -----------

###  current pred deviations for O. koreanus
o.cur.dev <- c(o.clim.dev, o.full.dev)
names(o.cur.dev) = c('Climate only', 'Full')

gplot(o.cur.dev) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 1, ncol = 2) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Deviation',
                       breaks = c(0.1, 0.8),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

ggsave('plots/Pred deviation/O.koreanus_current_pred_deviation.png', width = 15, height = 10, dpi = 800, units = 'cm')


###  current pred deviations for K. koreana
k.cur.dev <- c(k.clim.dev, k.full.dev)
names(k.cur.dev) = c('Climate only', 'Full')

gplot(k.cur.dev) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 1, ncol = 2) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Deviation',
                       breaks = c(0.1, 0.9),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

ggsave('plots/Pred deviation/K.koreana_current_pred_deviation.png', width = 15, height = 10, dpi = 800, units = 'cm')


###  hindcast pred deviations for O. koreanus
gplot(o.hind.dev) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, ncol = 5, nrow = 1) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Deviation',
                       breaks = c(0.1, 0.5),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

ggsave('plots/Pred deviation/O.koreanus_hindcast_deviation.png', width = 40, height = 15, dpi = 800, units = 'cm')


###  hindcast pred deviations for K. koreana
gplot(k.hind.dev) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, ncol = 5, nrow = 1) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Deviation',
                       breaks = c(0.1, 0.9),
                       labels = c('Low', 'High')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

ggsave('plots/Pred deviation/K.koreana_hindcast_deviation.png', width = 40, height = 15, dpi = 800, units = 'cm')

