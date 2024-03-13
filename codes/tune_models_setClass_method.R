####  model tuning function == if you want to use "setClass" instead of making a list.......
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