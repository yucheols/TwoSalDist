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
# make function
format_data <- function(data, clim_data) {
  b12 <- as.data.frame(data$b12$sim$D)
  b21 <- as.data.frame(data$b21$sim$D)
  
  b12$type = 'b12'
  b21$tpye = 'b21'
  
  colnames(b12) = c('D', 'Type')
  colnames(b21) = c('D', 'Type')
  
  results <- rbind(b12, b21)
  results$clim_data = clim_data
  return(results)
}

# formatting
w.full.dat <- format_data(data = w.full, clim_data = 'WorldClim full')  # WorldClim full
w.trim.dat <- format_data(data = w.trim, clim_data = 'WorldClim trimmed')  # WorldClim trimmed
c.full.dat <- format_data(data = c.full, clim_data = 'CHELSA full')     # CHELSA full
c.trim.dat <- format_data(data = c.trim, clim_data = 'CHELSA trimmed')     # CHELSA trimmed

# combine data
all_data <- rbind(w.full.dat, w.trim.dat, c.full.dat, c.trim.dat)
all_data$clim_data <- factor(all_data$clim_data, levels = c('WorldClim full', 'WorldClim trimmed', 'CHELSA full', 'CHELSA trimmed'))

# empirical D values
emp.d <- data.frame(var = c('WorldClim_full', 'WorldClim_trimmed', 'CHELSA_full', 'CHELSA_trimmed'),
                    val = c(w.full$b12$obs$D, w.trim$b12$obs$D, c.full$b12$obs$D, c.trim$b12$obs$D),
                    clim_data = c('WorldClim full', 'WorldClim trimmed', 'CHELSA full', 'CHELSA trimmed')) 

# plot
all_data %>%
  ggplot(aes(x = D, fill = Type, color = Type)) +
  facet_wrap(~ clim_data) +
  geom_histogram(bins = 10, alpha = 0.3, linewidth = 1.0) +
  geom_vline(data = emp.d, aes(xintercept = val), color = 'cornflowerblue', linewidth = 1.2, linetype = 'longdash') +
  scale_fill_manual(values = c('#69b3a2', '#9966ff')) +
  scale_color_manual(values = c('#69b3a2', '#9966ff')) +
  xlab("Schoener's D") + ylab('Count') +
  theme_bw() +
  theme(strip.text = element_text(size = 13))
