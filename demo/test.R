#!/usr/bin/env Rscript
list.of.packages <- c("parallel", "ArgumentCheck", "optparse")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
install.packages(file.path(getwd(), "evolike",fsep = .Platform$file.sep), repos = NULL, type="source")
library(evolike)
library(parallel)
library(optparse)


option_list = list(
  make_option(c("-p", "--prot"), type="character", default=NULL, 
              help="protein: h1, h3, n1 or n2", metavar="character"),
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
  
  if (!(model %in% c("weibull", "exponential"))){
    ArgumentCheck::addError(
      msg = "valid 'model' options: weibull, exponential",
      argcheck = Check
    )
  }
  
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  
  
  ## load mygroups 
  data("mygroups")
  print (" Yaaaaaaap! ")