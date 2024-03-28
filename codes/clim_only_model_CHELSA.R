# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(raster)
library(dplyr)
library(ENMeval)

##### Part 13 ::: fit climate-only model ---------------------------------------------------------------------------------------------
####  prep data 
## load occs
o.occs <- read.csv('data/occs/Onychodactylus_koreanus.csv') %>% select('long', 'lat')
k.occs <- read.csv('data/occs/Karsenia_koreana.csv') %>% select('long', 'lat')

head(o.occs)
head(k.occs)

## load bg
bg <- read.csv('data/bg/set1/bg1_10000.csv') %>% select('long', 'lat')
head(bg)

## load envs