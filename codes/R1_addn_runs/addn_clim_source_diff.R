#####  Visualize climate value differences at presence points between data sources

# clear working environment
rm(list = ls(all.names = T))
gc()


# load packages
library(dplyr)
library(raster)
library(ggplot2)
library(ggpubr)


### function for boxplot data formatting 
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

### load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv')

head(o.occs)
head(k.occs)

### load clim
envs_wc <- raster::stack(list.files(path = 'data/masked/WorldClim', pattern = '.bil$', full.names = T)) 
envs_wc <- raster::stack(subset(envs_wc, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))

envs_ch <- raster::stack(list.files(path = 'data/masked/CHELSA', pattern = '.bil$', full.names = T)) 
envs_ch <- raster::stack(subset(envs_ch, c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15')))


### WorldClim data formatting for boxplot
# O.koreanus
o.wc.dat <- boxdata(sp.name = 'O.koreanus', envs = envs_wc, pts = o.occs[, -1])
o.wc.dat$source = 'WorldClim'
head(o.wc.dat)

# K.koreana
k.wc.dat <- boxdata(sp.name = 'K.koreana', envs = envs_wc, pts = k.occs[, -1])
k.wc.dat$source = 'WorldClim'
head(k.wc.dat)

### CHELSA data formatting
# O.koreanus
o.ch.dat <- boxdata(sp.name = 'O.koreanus', envs = envs_ch, pts = o.occs[, -1])
o.ch.dat[1:374, 1] <- o.ch.dat[1:374, 1]/10  # div. bio1 and bio4 by 10 to match the unit
o.ch.dat$source = 'CHELSA'
head(o.ch.dat)

# K.koreana
k.ch.dat <- boxdata(sp.name = 'K.koreana', envs = envs_ch, pts = k.occs[, -1])
k.ch.dat[1:274, 1] <- k.ch.dat[1:274, 1]/10
k.ch.dat$source = 'CHELSA'
head(k.ch.dat)


### bind data per species for plotting
comb.dat <- rbind(o.wc.dat, o.ch.dat, k.wc.dat, k.ch.dat)
head(comb.dat)

### reorder plotting order
comb.dat$var <- factor(comb.dat$var, levels = c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15'))
comb.dat$Species <- factor(comb.dat$Species, levels = c('O.koreanus', 'K.koreana'))
comb.dat$source <- factor(comb.dat$source, levels = c('WorldClim', 'CHELSA'))

### recode var and species names
comb.dat$var <- dplyr::recode_factor(comb.dat$var,
                                     'bio1' = 'Bio1 (°C)', 'bio4' = 'Bio4', 'bio12' = 'Bio12 (mm)', 
                                     'bio13' = 'Bio13 (mm)', 'bio14' = 'Bio14 (mm)', 'bio15' = 'Bio15')


comb.dat$Species <- dplyr::recode_factor(comb.dat$Species, 'O.koreanus' = 'O. koreanus', 'K.koreana' = 'K. koreana')

### plot
comb.dat %>%
  ggplot(aes(x = Species, y = val, fill = source, color = source)) +
  geom_boxplot(linewidth = 2.0, alpha = 0.4, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(0.4), alpha = 0.4, size = 2.0) +
  facet_wrap(~ var, scale = 'free', nrow = 2, ncol = 3) +
  xlab('Species') + ylab('Value') +
  theme_bw() +
  theme(axis.title = element_text(size = 28, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 21),
        axis.text.x = element_text(face = 'italic'),
        strip.text = element_text(size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        legend.key.size = unit(9, 'mm'),
        legend.position = 'top') +
  stat_compare_means(method = 'wilcox.test',
                     label = 'p.signif',
                     tip.length = 0.01,
                     show.legend = F,
                     size = 7,
                     vjust = 0.5)

## save
ggsave('plots/compare_clim_source_revised.png', width = 35, height = 40, dpi = 800, units = 'cm')        


