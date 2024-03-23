#################  hindcasting
library(raster)

# climate only model
k.models2 <- ENMevaluate(taxon.name = 'K.koreana', occs = k.occs[, -1], envs = dropLayer(envs, c('forest', 'slope')), bg = bg1_10000[, -1], 
                         tune.args = list(fc = 'Q', rm = 1.0), algorithm = 'maxent.jar', doClamp = T, partitions = 'user', user.grp = k.folds[[2]])

# lig
lig <- raster::stack(list.files(path = 'hindcast_layers/LIG_Otto-Bliesner_30s', pattern = '.bil$', full.names = T))
