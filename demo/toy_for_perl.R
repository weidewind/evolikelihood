#!/usr/bin/env Rscript
list.of.packages <- c("parallel", "ArgumentCheck", "optparse")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
#install.packages(file.path(getwd(), "evolike",fsep = .Platform$file.sep), repos = NULL, type="source")
#library(evolike)
library(parallel)
library(optparse)


option_list = list(
  make_option(c("-p", "--prot"), type="character", default=NULL, 
              help="protein: h1, h3, n1 or n2", metavar="character"),
  make_option(c("-i", "--init_method"), type="character", default="cluster", 
              help="initialization method: clusterization of parameters (cluster) or radomly chosen parameters (random) [default= %default]", metavar="character"),
  make_option(c("-t", "--trials"), type="integer", default=1, 
              help="number of em trials (ignored if number of categories is 1 or if initialization method is cluster)", metavar="integer"),
  make_option(c("-c", "--categories"), type="integer", default=3, 
              help="number of categories", metavar="integer"),
  make_option(c("-m", "--model"), type="character", default="weibull", 
              help="model distr: weibull or exponential", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

with (opt, {
  
  Check <- ArgumentCheck::newArgCheck()
  
  if (!(prot %in% c("h1", "h3", "n1", "n2"))){ 
    ArgumentCheck::addError(
      msg = "valid 'prot' options: h1, h3, n1, n2",
      argcheck = Check
    )
  }
  if (!(init_method %in% c("cluster", "random"))){
    ArgumentCheck::addError(
      msg = "valid 'init_method' options: cluster, random",
      argcheck = Check
    )
  }
  
  if (categories < 1){
    ArgumentCheck::addError(
      msg = "'categories' must be >= 1",
      argcheck = Check
    )
  }
  if (trials < 1){
    ArgumentCheck::addError(
      msg = "'trials' must be >= 1",
      argcheck = Check
    )
  }
  if (!(model %in% c("weibull", "exponential"))){
    ArgumentCheck::addError(
      msg = "valid 'model' options: weibull, exponential",
      argcheck = Check
    )
  }
  if (categories == 1 && trials > 1){
    trials = 1
    ArgumentCheck::addWarning(
      msg = "'categories' equals to 1, will make only one trial (trials option ignored)",
      argcheck = Check
    )
  }
  if (init_method == "cluster"  && trials > 1){
    trials = 1
    ArgumentCheck::addWarning(
      msg = "initialization by clustering, will make only one trial (trials option ignored)",
      argcheck = Check
    )
  }
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  print (getwd())
  #print (file.path(getwd(), "output","toys", model, paste(c(prot, "_", init_method, "_", categories, "_", trial), collapse=""),fsep = .Platform$file.sep))
  #sink (file.path(getwd(), "output","toys", model, paste(c(prot, "_", init_method, "_", categories, "_", trial), collapse=""),fsep = .Platform$file.sep))
  #print ("So far so good")
  #sink()
  
 