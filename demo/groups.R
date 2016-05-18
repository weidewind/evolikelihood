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

prot_groups <-mygroups[grep(paste(c("^", prot, "_"), collapse = ""), names(mygroups))]
categories <- 1

prot_data <-  read.csv(file.path(getwd(), "data", paste(c(prot,"_for_LRT.csv"), collapse=""),fsep = .Platform$file.sep),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)

count_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(count_cores)
clusterExport(cl, list("prot", "prot_groups", "splitted", "model"), envir = environment())
clusterCall(cl, function() library(evolike))

prot_groups_names <- names(prot_groups)
em_results_list <- parLapply(cl, prot_groups_names, function(group_name){
  sink (file.path(getwd(), "output","group_likelihood", model, prot, paste(c("logs_", group_name), collapse= ""), fsep = .Platform$file.sep))
  
  print (group_name)
  group <- prot_groups[group_name][[1]]
  group_pattern <- sapply(group, function (e){
    paste(c("^", e, "\\."), collapse = "")
  })
  grepper <- paste(group_pattern,collapse="|")
  group_nodes <- grep(grepper, names(splitted))
  group_splitted <- splitted[group_nodes]
  
  compl <- !(seq(1,length(splitted),1) %in% group_nodes)
  compl_splitted <- splitted[compl]
  
  params <-parameters(group_splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
  compl_params <-parameters(compl_splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
  
  sink() 
  
  sink (file.path(getwd(), "output","group_likelihood", model, prot, group_name, fsep = .Platform$file.sep))
  em_results <- em_procedure(data=group_splitted, params=params, model = model,  cluster.number= 1, init_method = "cluster", mutation_position = "middle",  filtering = "single", trace = FALSE)
  sink() 
  sink (file.path(getwd(), "output","group_likelihood", model, prot, paste(c(group_name, "_complement"), collapse= ""), fsep = .Platform$file.sep))
  em_results <- em_procedure(data=compl_splitted, params=compl_params, model = model,  cluster.number= 1, init_method = "cluster", mutation_position = "middle",  filtering = "single", trace = FALSE)
  sink() 
  
  
  em_results
})




stopCluster(cl)