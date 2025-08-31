############   plot CHELSA-optimized prediction results

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(terra)
library(rgdal) 
library(rasterVis)
library(ggplot2)
library(patchwork)

# load korea polygon
kor <- readOGR('data/polygons/kor_mer.shp')


#####  plot O.koreanus CHELSA current maps

###  plot continuous
# continuous clim only
o.clim.ch <- rast('tuning_experiments/preds_addn/O.koreanus/CHELSA_clim_only/preds/bg1_10000.tif')
plot(o.clim.ch)
names(o.clim.ch) = 'Continuous'

# plot
o.cont.ch.plot <- gplot(o.clim.ch) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))
  
###  plot binary
# binary clim only
o.clim.ch.bin <- rast('tuning_experiments/preds_addn/O.koreanus/bin_both/O.koreanus_addn_clim_only_bin.tif')
plot(o.clim.ch.bin)
names(o.clim.ch.bin) = 'Binary'

# plot
o.bin.ch.plot <- gplot(o.clim.ch.bin) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = rev(terrain.colors(1000)),
                       na.value = 'NA') +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.position = 'none',
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

### arrange plots
(o.cont.ch.plot | o.bin.ch.plot)

### export
ggsave('plots/CHELSA tuning/O.koreanus_CHELSA_tune.png', width = 20, height = 10, dpi = 800, units = 'cm')



#####  plot K.koreana CHELSA current maps

###  plot continuous
# continuous clim only
k.clim.ch <- rast('tuning_experiments/preds_addn/K.koreana/CHELSA_clim_only/preds/bg1_10000.tif')
plot(k.clim.ch)
names(k.clim.ch) = 'Continuous'

# plot
k.cont.ch.plot <- gplot(k.clim.ch) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))


###  plot binary
# binary clim only
k.clim.ch.bin <- rast('tuning_experiments/preds_addn/K.koreana/bin_both/K.koreana_addn_clim_only_bin.tif')
plot(k.clim.ch.bin)
names(k.clim.ch.bin) = 'Binary'

# plot
k.bin.ch.plot <- gplot(k.clim.ch.bin) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = rev(terrain.colors(1000)),
                       na.value = NA) + 
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.position = 'none',
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 12))

### arrange plots
(k.cont.ch.plot | k.bin.ch.plot)

### export
ggsave('plots/CHELSA tuning/K.koreana_CHELSA_tune.png', width = 20, height = 10, dpi = 800, units = 'cm')


#####  plot O.koreanus CHELSA hindcast maps
# load prediction rasters
o.hinds <- rast(list.files(path = 'tuning_experiments/preds_addn/O.koreanus/CHELSA_clim_only/hindcast/', pattern = '.tif', full.names = T))
o.hinds <- o.hinds[[c(5, 4, 2, 1, 3)]] 
print(o.hinds)

# plot
gplot(o.hinds) + 
  geom_tile(aes(fill = value)) +
  coord_equal(expand = F) +
  facet_wrap(~ variable, nrow = 1, ncol = 5) +
  scale_x_continuous(breaks = c(122, 128, 134)) +
  scale_fill_gradientn(colors = c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 17),
        legend.title = element_text(size = 17, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 17, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 15))

# export
ggsave('plots/CHELSA tuning/O.koreanus_hindcast_CHELSA_tune.png', width = 40, height = 7, dpi = 800, units = 'cm')


#####  plot K.koreana CHELSA hindcast maps
# load prediction rasters
k.hinds <- rast(list.files(path = 'tuning_experiments/preds_addn/K.koreana/CHELSA_clim_only/hindcast/', pattern = '.tif', full.names = T))
k.hinds <- k.hinds[[c(5,4,2,1,3)]]
print(k.hinds)

# plot
gplot(k.hinds) +
  geom_tile(aes(fill = value)) +
  coord_equal(expand = F) +
  facet_wrap(~ variable, nrow = 1, ncol = 5) +
  scale_x_continuous(breaks = c(122, 128, 134)) +
  scale_fill_gradientn(colors = c('#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#4f05d7'),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 17),
        legend.title = element_text(size = 17, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 17, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 15))

# export
ggsave('plots/CHELSA tuning/K.koreana_hindcast_CHELSA_tune.png', width = 40, height = 7, dpi = 800, units = 'cm')

