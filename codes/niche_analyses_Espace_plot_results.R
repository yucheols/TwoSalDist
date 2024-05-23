#### clean up working env
rm(list = ls(all.names = T))
gc()

#### plot results
library(ggplot2)
library(dplyr)

##### Step 4 ::: plot full run results -------------------------------------------------------------------------------------------------
### load result
w.full <- readRDS('niche_analyses/WorldClim/E-space/output/full_run_mx.rds')
w.trim <- readRDS('niche_analyses/WorldClim/E-space/output/trimmed_run_mx.rds')
c.full <- readRDS('niche_analyses/CHELSA/E-space/output/full_run_mx.rds')
c.trim <- readRDS('niche_analyses/CHELSA/E-space/output/trimmed_run_mx.rds')

### data formatting
# WorldClim full
full.b12.d <- as.data.frame(full.run.mx$b12$sim$D)
full.b21.d <- as.data.frame(full.run.mx$b21$sim$D)

full.b12.d$type = 'b12'
full.b21.d$type = 'b21'

colnames(full.b12.d) = c('D', 'type')
colnames(full.b21.d) = c('D', 'type')

full.results <- rbind(full.b12.d, full.b21.d)
head(full.results)

# plot
full.results %>%
  ggplot(aes(x = D, fill = type)) +
  geom_histogram()
