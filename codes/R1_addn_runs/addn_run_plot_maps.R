############   plot prediction results

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(terra)
library(rgdal) 
library(rasterVis)
library(ggplot2)

# load korea polygon
kor <- readOGR('data/polygons/kor_mer.shp')

#####  plot O.koreanus CHELSA current maps
o.cur.ch <- rast('tuning_experiments/preds_addn/O.koreanus/CHELSA_clim_only/preds/bg1_10000.tif')



#####  plot O.koreanus CHELSA hindcast maps
# load prediction rasters
o.hinds <- rast(list.files(path = 'tuning_experiments/preds_addn/O.koreanus/CHELSA_clim_only/hindcast/', pattern = '.tif', full.names = T))
o.hinds <- o.hinds[[c(5, 4, 2, 1, 3)]] 
print(o.hinds)

# plot
gplot(o.hinds) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 1, ncol = 5) +
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

# export
ggsave('plots/hindcast/O.koreanus_hindcast_CHELSA_tune.png', width = 40, height = 7, dpi = 800, units = 'cm')


#####  plot K.koreana CHELSA hindcast maps
# load prediction rasters
k.hinds <- rast(list.files(path = 'tuning_experiments/preds_addn/K.koreana/CHELSA_clim_only/hindcast/', pattern = '.tif', full.names = T))
k.hinds <- k.hinds[[c(5,4,2,1,3)]]
print(k.hinds)

# plot
gplot(k.hinds) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, nrow = 1, ncol = 5) +
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

# export
ggsave('plots/hindcast/K.koreana_hindcast_CHELSA_tune.png', width = 40, height = 7, dpi = 800, units = 'cm')

