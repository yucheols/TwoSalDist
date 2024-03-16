##### Niche analyses in E-space using humboldt
# setwd
setwd('niche_analyses')

# clean up working env
rm(list = ls(all.names = T))
gc()

# load packages
library(humboldt)
library(raster)

##### Step 1 ::: prep env data ----------------------------------------------------------------------------------------------------
# use WorldClim 1km layers
# use the selected subset of variables used for modeling
e.data <- raster::stack(list.files(path = 'data/WorldClim', pattern = '.bil$', full.names = T))
e.data <- raster::stack(subset(e.data, c('bio1', 'bio2', 'bio3', 'bio4', 'bio12', 'bio13', 'bio14', 'forest', 'slope')))
plot(e.data[[1]])

####  Convert input rasters to 'humboldt' format for use as an environment file
e.points <- rasterToPoints(e.data[[1]], fun = NULL, spatial = F)
print(e.points)

## subset only the x and y data
e.points <- e.points[, c(1,2)]
head(e.points)

## Extract values to points from rasters
ras.val <- data.frame(raster::extract(e.data, e.points))
head(ras.val)

## merge sampled data to input
env.data <- cbind(e.points, ras.val)
head(env.data)

## remove NAs and make sure all variables are imported as numbers
env.data <- humboldt.scrub.env(env.data)

## save the file as '.csv' for future analyses 
write.csv(env.data, 'data/niche_analyses/env_data_wc_subset.csv')


#####  Step 2 ::: load occs data and format it for analysis in humboldt  ---------------------------------------------------------------------------------
#####  colnames must be ( 'sp', 'x', 'y' )

## occs.sp1 :: O.koreanus
occs.sp1 <- read.csv('data/occs/Onychodactylus_koreanus.csv')
colnames(occs.sp1) = c('sp', 'x', 'y')
head(occs.sp1)


## occs.sp2 :: K.koreana
occs.sp2 <- read.csv('data/occs/Karsenia_koreana.csv')
colnames(occs.sp2) = c('sp', 'x', 'y')
head(occs.sp2)


#####  Step 3 ::: Niche Overlap and Niche Divergence Tests -------------------------------------------------------------------------------------------------

## run it first with full environmental for background tests and equivalence statistic 
## (total equivalence or divergence in current distributions)  
full.run.mx <- humboldt.doitall(inname = 'full_run_mx', env1 = env.data, env2 = env.data, sp1 = occs.sp1, sp2 = occs.sp2,
                                rarefy.dist = 1, rarefy.units = 'km', env.reso = 0.0083333333, reduce.env = 0,
                                non.analogous.environments = 'YES', correct.env = T, env.trim = T, env.trim.type = 'RADIUS',
                                trim.buffer.sp1 = 50, trim.buffer.sp2 = 50, color.ramp = 5,
                                pcx = 1, pcy = 2, col.env = e.var, e.var = c(3:ncol(env.data)),
                                R = 100, kern.smooth = 1, e.reps = 1000, b.reps = 1000, nae = 'YES', thresh.espace.z = 0.001,
                                p.overlap = T, p.boxplot = T, p.scatter = T, run.silent = F, ncores = 'all')


## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche divergence. 
## For comparing results, change only the following model parameters: reduce.env, non.analogous.environmental, env.trim, nae
trimmed.run.mx <- humboldt.doitall(inname = 'trimmed_run_mx', env1 = env.data, env2 = env.data, sp1 = occs.sp1, sp2 = occs.sp2,
                                   rarefy.dist = 1, rarefy.units = 'km', env.reso = 0.0083333333, reduce.env = 2, reductype = 'PCA',
                                   non.analogous.environments = 'NO', correct.env = T, env.trim = T, env.trim.type = 'RADIUS',
                                   trim.buffer.sp1 = 50, trim.buffer.sp2 = 50, color.ramp = 5,
                                   pcx = 1, pcy = 2, col.env = e.var, e.var = c(3:ncol(env.data)),
                                   R = 100, kern.smooth = 1, e.reps = 1000, b.reps = 1000, nae = 'NO', thresh.espace.z = 0.001,
                                   p.overlap = T, p.boxplot = T, p.scatter = T, run.silent = F, ncores = 'all')


## Save outputs
#dir.create('E-space/output')

#saveRDS(full.run.mx, 'E-space/output/full_run_mx.rds')
#saveRDS(trimmed.run.mx, 'E-space/output/trimmed_run_mx.rds')
