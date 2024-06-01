#####  Visualize climate value differences at presence points between data sources

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(dplyr)
library(raster)
library(ggplot2)

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
  ggplot(aes(x = source, y = val, fill = Species, color = Species)) +
  geom_boxplot(linewidth = 1.0, alpha = 0.4, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(0.4), alpha = 0.4) +
  facet_wrap(~ var, scale = 'free', nrow = 2, ncol = 3) +
  scale_fill_manual(values = c('#6495ED', '#ffe600')) +
  scale_color_manual(values = c('#6495ED', '#ffe600')) +
  xlab('Source') + ylab('Value') +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14, face = 'italic'),
        legend.position = 'top')

## save
ggsave('plots/compare_clim_source.png', width = 20, height = 25, dpi = 800, units = 'cm')        


##### run some stats analysis

# turn off scientific notation 
options(scipen = 999)

#### O.koreanus
## filter by variable
# WorldClim
o.wc.b1 <- o.wc.dat %>% filter(var == 'bio1')
o.wc.b4 <- o.wc.dat %>% filter(var == 'bio4')
o.wc.b12 <- o.wc.dat %>% filter(var == 'bio12')
o.wc.b13 <- o.wc.dat %>% filter(var == 'bio13')
o.wc.b14 <- o.wc.dat %>% filter(var == 'bio14')
o.wc.b15 <- o.wc.dat %>% filter(var == 'bio15')

# CHELSA
o.ch.b1 <- o.ch.dat %>% filter(var == 'bio1')
o.ch.b4 <- o.ch.dat %>% filter(var == 'bio4')
o.ch.b12 <- o.ch.dat %>% filter(var == 'bio12')
o.ch.b13 <- o.ch.dat %>% filter(var == 'bio13')
o.ch.b14 <- o.ch.dat %>% filter(var == 'bio14')
o.ch.b15 <- o.ch.dat %>% filter(var == 'bio15')

## run tests
o.mann.b1 <- wilcox.test(x = o.wc.b1$val, y = o.ch.b1$val, alternative = 'two.sided')          # bio1
o.mann.b4 <- wilcox.test(x = o.wc.b4$val, y = o.ch.b4$val, alternative = 'two.sided')          # bio4
o.mann.b12 <- wilcox.test(x = o.wc.b12$val, y = o.ch.b12$val, alternative = 'two.sided')       # bio12
o.mann.b13 <- wilcox.test(x = o.wc.b13$val, y = o.ch.b13$val, alternative = 'two.sided')       # bio13
o.mann.b14 <- wilcox.test(x = o.wc.b14$val, y = o.ch.b14$val, alternative = 'two.sided')       # bio14
o.mann.b15 <- wilcox.test(x = o.wc.b15$val, y = o.ch.b15$val, alternative = 'two.sided')       # bio15

## make a dataframe of results
o.mann <- data.frame(W = c(o.mann.b1$statistic, o.mann.b4$statistic, o.mann.b12$statistic, o.mann.b13$statistic, o.mann.b14$statistic, o.mann.b15$statistic), 
                     p = c(o.mann.b1$p.value, o.mann.b4$p.value, o.mann.b12$p.value, o.mann.b13$p.value, o.mann.b14$p.value, o.mann.b15$p.value),
                     variable = c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15'),
                     species = 'O.koreanus')

print(o.mann)

## save results
write.csv(o.mann, 'clim_source_compare/O.koreanus_MannWhitneyU.csv')


#### K.koreana
## filter by variable
# WorldClim
k.wc.b1 <- k.wc.dat %>% filter(var == 'bio1')
k.wc.b4 <- k.wc.dat %>% filter(var == 'bio4')
k.wc.b12 <- k.wc.dat %>% filter(var == 'bio12')
k.wc.b13 <- k.wc.dat %>% filter(var == 'bio13')
k.wc.b14 <- k.wc.dat %>% filter(var == 'bio14')
k.wc.b15 <- k.wc.dat %>% filter(var == 'bio15')

# CHELSA
k.ch.b1 <- k.ch.dat %>% filter(var == 'bio1')
k.ch.b4 <- k.ch.dat %>% filter(var == 'bio4')
k.ch.b12 <- k.ch.dat %>% filter(var == 'bio12')
k.ch.b13 <- k.ch.dat %>% filter(var == 'bio13')
k.ch.b14 <- k.ch.dat %>% filter(var == 'bio14')
k.ch.b15 <- k.ch.dat %>% filter(var == 'bio15')

## run tests
k.mann.b1 <- wilcox.test(x = k.wc.b1$val, y = k.ch.b1$val, alternative = 'two.sided')          # bio1
k.mann.b4 <- wilcox.test(x = k.wc.b4$val, y = k.ch.b4$val, alternative = 'two.sided')          # bio4
k.mann.b12 <- wilcox.test(x = k.wc.b12$val, y = k.ch.b12$val, alternative = 'two.sided')       # bio12
k.mann.b13 <- wilcox.test(x = k.wc.b13$val, y = k.ch.b13$val, alternative = 'two.sided')       # bio13
k.mann.b14 <- wilcox.test(x = k.wc.b14$val, y = k.ch.b14$val, alternative = 'two.sided')       # bio14
k.mann.b15 <- wilcox.test(x = k.wc.b15$val, y = k.ch.b15$val, alternative = 'two.sided')       # bio15

## make a dataframe of results
k.mann <- data.frame(W = c(k.mann.b1$statistic, k.mann.b4$statistic, k.mann.b12$statistic, k.mann.b13$statistic, k.mann.b14$statistic, k.mann.b15$statistic), 
                     p = c(k.mann.b1$p.value, k.mann.b4$p.value, k.mann.b12$p.value, k.mann.b13$p.value, k.mann.b14$p.value, k.mann.b15$p.value),
                     variable = c('bio1', 'bio4', 'bio12', 'bio13', 'bio14', 'bio15'),
                     species = 'K.koreana')

print(k.mann)

## save results
write.csv(k.mann, 'clim_source_compare/K.koreana_MannWhitneyU.csv')
