#### clean up working env
rm(list = ls(all.names = T))
gc()

#### plot results
library(ggplot2)
library(dplyr)
library(ggpubr)

##### Step 4 ::: plot full run results -------------------------------------------------------------------------------------------------
### load result
w.full <- readRDS('niche_analyses/WorldClim/E-space/output/full_run_mx.rds')
w.trim <- readRDS('niche_analyses/WorldClim/E-space/output/trimmed_run_mx.rds')
c.full <- readRDS('niche_analyses/CHELSA/E-space/output/full_run_mx.rds')
c.trim <- readRDS('niche_analyses/CHELSA/E-space/output/trimmed_run_mx.rds')

### data formatting
# make function
format_espace_data <- function(data, clim_data) {
  sim <- as.data.frame(data$eqiv$sim$D)
  b12 <- as.data.frame(data$b12$sim$D)
  b21 <- as.data.frame(data$b21$sim$D)
  
  sim$type = 'Identity'
  b12$type = 'b12'
  b21$tpye = 'b21'
  
  colnames(sim) = c('D', 'Type')
  colnames(b12) = c('D', 'Type')
  colnames(b21) = c('D', 'Type')
  
  id_test_results <- sim
  id_test_results$clim_data = clim_data
  
  bg_test_results <- rbind(b12, b21)
  bg_test_results$clim_data = clim_data
  
  results <- list(id_test_results, bg_test_results)
  
  return(results)
}

# formatting
w.full.dat <- format_espace_data(data = w.full, clim_data = 'WorldClim full')  # WorldClim full
w.trim.dat <- format_espace_data(data = w.trim, clim_data = 'WorldClim trimmed')  # WorldClim trimmed
c.full.dat <- format_espace_data(data = c.full, clim_data = 'CHELSA full')     # CHELSA full
c.trim.dat <- format_espace_data(data = c.trim, clim_data = 'CHELSA trimmed')     # CHELSA trimmed

#### identity test 
# combine data
id.test.result <- rbind(w.full.dat[[1]], w.trim.dat[[1]], c.full.dat[[1]], c.trim.dat[[1]])

# empirical D values
id.emp.d <- data.frame(var = c('WorldClim_full', 'WorldClim_trimmed', 'CHELSA_full', 'CHELSA_trimmed'),
                       val = c(w.full$eqiv$obs$D, w.trim$eqiv$obs$D, c.full$eqiv$obs$D, c.trim$eqiv$obs$D),
                       clim_data = c('WorldClim full', 'WorldClim trimmed', 'CHELSA full', 'CHELSA trimmed'))

# plot
id.plot <- id.test.result %>% 
  ggplot(aes(x = D, fill = Type, color = Type)) + 
  geom_histogram(bins = 10, alpha = 0.3, linewidth = 1.0) +
  facet_wrap(~ factor(clim_data, levels = c('WorldClim full', 'WorldClim trimmed', 'CHELSA full', 'CHELSA trimmed')), nrow = 1, ncol = 4) + 
  geom_vline(data = id.emp.d, aes(xintercept = val), color = 'black', linewidth = 1.0, linetype = 'longdash') +
  scale_fill_manual(values = 'cornflowerblue') +
  scale_color_manual(values = 'cornflowerblue') +
  xlab("Schoener's D") + ylab('Count') +
  theme_bw() +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20))) 


#### background test
# combine data
bg.test.data <- rbind(w.full.dat[[2]], w.trim.dat[[2]], c.full.dat[[2]], c.trim.dat[[2]])

# recode the "Type" column
bg.test.data$Type <- recode_factor(bg.test.data$Type,
                                   'b12' = 'Sp2 vs Sp1 background',
                                   'b21' = 'Sp1 vs Sp2 background')

# empirical D values
bg.emp.d <- data.frame(var = c('WorldClim_full', 'WorldClim_trimmed', 'CHELSA_full', 'CHELSA_trimmed'),
                       val = c(w.full$b12$obs$D, w.trim$b12$obs$D, c.full$b12$obs$D, c.trim$b12$obs$D),
                       clim_data = c('WorldClim full', 'WorldClim trimmed', 'CHELSA full', 'CHELSA trimmed')) 

# plot
bg.plot <- bg.test.data %>%
  ggplot(aes(x = D, fill = Type, color = Type)) +
  geom_histogram(bins = 10, alpha = 0.3, linewidth = 1.0) +
  facet_wrap(~ factor(clim_data, levels = c('WorldClim full', 'WorldClim trimmed', 'CHELSA full', 'CHELSA trimmed')), nrow = 1, ncol = 4) +
  geom_vline(data = bg.emp.d, aes(xintercept = val), color = 'black', linewidth = 1.0, linetype = 'longdash') +
  scale_fill_manual(values = c('#69b3a2', '#9966ff')) +
  scale_color_manual(values = c('#69b3a2', '#9966ff')) +
  xlab("Schoener's D") + ylab('Count') +
  theme_bw() +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20))) +
  guides(fill = guide_legend(reverse = T),
         color = guide_legend(reverse = T))


### combine plots
espace_fig <- ggarrange(id.plot, bg.plot,
                        font.label = list(size = 16),
                        labels = c('A', 'B'),
                        ncol = 1, nrow = 2)

### save
ggsave('plots/niche_analyses_Espace.png', width = 25, height = 20, dpi = 800, units = 'cm')  
