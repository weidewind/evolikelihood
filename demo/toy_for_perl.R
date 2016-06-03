#!/usr/bin/env Rscript
list.of.packages <- c("parallel", "ArgumentCheck", "optparse")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
install.packages(file.path(getwd(), "evolike",fsep = .Platform$file.sep), repos = NULL, type="source")
#wd <- getwd()
#setwd("..")
#parent <- getwd()
#setwd(wd)
#install.packages(file.path(parent, "evolike", fsep = .Platform$file.sep), repos = NULL, type="source")
library(evolike)
library(parallel)
library(optparse)


option_list = list(
  make_option(c("-p", "--prot"), type="character", default=NULL, 
              help="protein: h1, h3, n1 or n2", metavar="character"),
  make_option(c("-i", "--init_method"), type="character", default="cluster", 
              help="initialization method: clusterization of parameters (cluster) or radomly chosen parameters (random) [default= %default]", metavar="character"),
  make_option(c("-t", "--trial"), type="integer", default=1, 
              help="pid of trial", metavar="integer"),
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

  if (!(model %in% c("weibull", "exponential"))){
    ArgumentCheck::addError(
      msg = "valid 'model' options: weibull, exponential",
      argcheck = Check
    )
  }

  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  #print (getwd())
  #print (file.path(getwd(), "output","toys", model, paste(c(prot, "_", init_method, "_", categories, "_", trial), collapse=""),fsep = .Platform$file.sep))
  #sink (file.path(getwd(), "output","toys", model, paste(c(prot, "_", init_method, "_", categories, "_", trial), collapse=""),fsep = .Platform$file.sep))
  #print ("So far so good")
  #sink()
  print (getwd())
  print (em_procedure)
  prot_data <-  read.csv(file.path(getwd(), "data", paste(c(prot,"_for_LRT.csv"), collapse=""),fsep = .Platform$file.sep),stringsAsFactors=FALSE)  
  splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
  params <-parameters(splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
  
  sink (file.path(getwd(), "output","toys", model, prot, paste(c(prot, "_", init_method, "_", categories, "_", trial), collapse=""),fsep = .Platform$file.sep))
  trackfile <-file.path(getwd(), "output","toys", model, prot, paste(c(prot, "_", init_method, "_", categories, "_", trial, "_track_"), collapse=""),fsep = .Platform$file.sep)
  em_results <- em_procedure(data=splitted, params=params, model = model, iter = 1000, cluster.number= categories, init_method = init_method, mutation_position = "middle",  filtering = "single", trace = FALSE, trackfile = trackfile, trackcount = 1)
  sink() 
  
})
  
 