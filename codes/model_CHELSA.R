##### Data prep for current ENM with CHELSA data 
# set seed 
set.seed(321)

# load packages
library(ENMeval)
library(dplyr)

# check data
print(envs)

head(o.occs)
head(k.occs)

head(bg1_5000)
head(bg2_5000)


#####  Part 6 ::: model testing  ---------------------------------------------------------------------------------------------
# set output class
setClass(Class = 'ENMtuning',
         representation(output = 'list',
                        models = 'list',
                        preds = 'list'))

# automate model tuning 
test_models <- function(taxon.name, occs, envs, bg.list, tune.args, partitions, user.grp) {
  output <- list()
  models <- list()
  preds <- list()
  
  for (i in 1:length(bg.list)) {
    
    # make models
    eval <- ENMeval::ENMevaluate(taxon.name = taxon.name, occs = occs, envs = envs,
                                 bg = bg.list[[i]], tune.args = tune.args, partitions = partitions,
                                 user.grp = user.grp[[i]], doClamp = T, algorithm = 'maxent.jar', parallel = T,
                                 parallelType = 'doSNOW')
    
    # get results
    eval.res <- ENMeval::eval.results(eval)
    
    # get optimal parameter combinations
    opt.param <- eval.res %>% dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
      dplyr::filter(auc.diff.avg == min(auc.diff.avg)) %>%
      dplyr::filter(auc.val.avg == max(auc.val.avg))
    
    output[[i]] <- opt.param
    
    # get optimal model per iteration
    opt.model <- ENMeval::eval.models(eval)[[opt.param$tune.args]]
    models[[i]] <- opt.model
    
    # get optimal predictions per iteration
    opt.pred <- ENMeval::eval.predictions(eval)[[opt.param$tune.args]]
    preds[[i]] <- opt.pred
  }
  return(new('ENMtuning', output = output, models = models, preds = preds))
}


### prep inputs
bg.list <- list(bg1_5000[, c('long', 'lat')], bg1_10000[, c('long', 'lat')], bg1_15000[, c('long', 'lat')],
                bg2_5000[, c('long', 'lat')], bg2_10000[, c('long', 'lat')], bg2_15000[, c('long', 'lat')])

tune.args <- list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'), 
                  rm = seq(0.5,5, by = 0.5))


### O. koreanus model testing run
# run
o.models <- test_models(taxon.name = 'O.koreanus', occs = o.occs[, c(2,3)], envs = envs, 
                        bg.list = bg.list, tune.args = tune.args, partitions = c('user'), user.grp = o.folds)

# look at results


# look at predictions


### K. koreana model testing run
k.models <- test_models(taxon.name = 'K.koreana', occs = k.occs[, c(2,3)], envs = envs, 
                        bg.list = bg.list, tune.args = tune.args, partitions = c('user'), user.grp = k.folds)

# look at results


# look at predictions
