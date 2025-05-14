#####  sptially visulize WorldClim vs. CHELSA climate layer differences
#####  use raw difference to visualize direction and magnitude of difference // visualize systematic bias or directional differences between datasets

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(terra)
library(rgdal) 
library(rasterVis)
library(ggplot2)
library(ggpubr)
library(pals)

# load korea polygon
kor <- readOGR('data/polygons/kor_mer.shp')


#####  load WorldClim layers
wc <- rast(list.files(path = 'data/masked/WorldClim/', pattern = '.bil$', full.names = T))
wc <- terra::subset(wc, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15'))
print(wc)

#####  load CHELSA layers
ch <- rast(list.files(path = 'data/masked/CHELSA/', pattern = '.bil$', full.names = T))
ch <- terra::subset(ch, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15'))
print(ch)

#####  compare units == need to divide CHELSA bio1 and bio4 by 10
ENMwrap::unit_check(ref.env = raster::stack(wc), proj.env = raster::stack(ch), n = 10000)

ch[['bio1']] <- ch[['bio1']]/10
ch[['bio4']] <- ch[['bio4']]/10

plot(ch[['bio1']])
plot(ch[['bio4']])

#####  calculate deviation
clim.dev <- wc - ch
plot(clim.dev)


#####  plot
clim.dev <- clim.dev[[c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')]]
names(clim.dev) = c('Bio 1 (°C)', 'Bio 4', 'Bio 12 (mm)', 'Bio 13 (mm)', 'Bio 14 (mm)', 'Bio 15')
print(clim.dev)

# bio1 plot
bio1_plot <- gplot(clim.dev[[1]]) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Difference') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text = element_text(size = 12))


# bio4 plot
bio4_plot <- gplot(clim.dev[[2]]) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Difference') +
  #xlab('Longitude (°)') + ylab('Latitude (°)') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text = element_text(size = 12))


# bio12 plot
bio12_plot <- gplot(clim.dev[[3]]) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Difference') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text = element_text(size = 12))


# bio13 plot
bio13_plot <- gplot(clim.dev[[4]]) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Difference') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text = element_text(size = 12))


# bio14
bio14_plot <- gplot(clim.dev[[5]]) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Difference') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text = element_text(size = 12))


# bio15
bio15_plot <- gplot(clim.dev[[6]]) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors =  rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Difference') +
  geom_polygon(data = kor, aes(x = long, y = lat, group = group), color = 'black', linewidth = 0.5, linetype = 'solid', fill = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, margin = margin(b = 10)),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text = element_text(size = 12))

###  combine plots 
comb_clim_plot <- ggarrange(bio1_plot, bio4_plot, bio12_plot, bio13_plot, bio14_plot, bio15_plot, ncol = 3, nrow = 2, align = 'hv')
annotate_figure(comb_clim_plot, 
                left = text_grob('Latitude (°)', rot = 90, size = 14, face = 'bold'),
                bottom = text_grob('Longitude (°)', size = 14, face = 'bold'))

# export
ggsave('plots/Clim deviation/clim_layers_deviation.png', width = 30, height = 15, dpi = 800, units = 'cm')




###  scale the deviation to 0-1 range
# automate
rast_scaler <- function(input_rast_1 , input_rast_2) {
  diff_rast <- rast()
  
  for (i in 1:nlyr(input_rast_1)) {
    diff <- input_rast_1[[i]] - input_rast_2[[i]]
    rast_min <- global(diff, min, na.rm = T)[1,1]
    rast_max <- global(diff, max, na.rm = T)[1,1]
    diff_scl <- (diff - rast_min) / (rast_max - rast_min)
    diff_rast <- c(diff_rast, diff_scl)
  }
  names(diff_rast) = names(input_rast_1)
  return(diff_rast)
}

# run
scale_clim <- rast_scaler(input_rast_1 = wc, input_rast_2 = ch)
plot(scale_clim)


# plot scaled layers
names(scale_clim) = c('Bio 1 (°C)', 'Bio 4', 'Bio 12 (mm)', 'Bio 13 (mm)', 'Bio 14 (mm)', 'Bio 15')

gplot(scale_clim) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colors = rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Difference',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
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

# export
ggsave('plots/Clim deviation/clim_layers_deviation_scaled.png', width = 30, height = 15, dpi = 800, units = 'cm')

