##########  plot hindcasting prediction maps
# clean up working env
rm(list = ls(all.names = T))
gc()

# load packages
library(rasterVis)
library(terra)
library(ggplot2)

### import polygon
poly <- rgdal::readOGR('data/polygons/kor_mer.shp')

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
