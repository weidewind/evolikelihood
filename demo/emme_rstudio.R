#!/usr/bin/env Rscript
list.of.packages <- c("parallel", "ArgumentCheck", "optparse")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
install.packages(file.path(getwd(), "evolike",fsep = .Platform$file.sep), repos = NULL, type="source")
library(evolike)
library(parallel)
library(optparse)


prot <- "n2"
model <- "weibull"
trials <- 20
init_method <- "random"
categories <- 2
diff_cats <-FALSE # parallel execution of cluster method for several categories, if true. parallel execution for several trials,if false

prot_data <-  read.csv(file.path(getwd(), "data", paste(c(prot,"_for_LRT.csv"), collapse=""),fsep = .Platform$file.sep),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
params <-parameters(splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)


  
if (diff_cats){
  cats <- seq(1, categories, 1)
}

# Calculate the number of cores
count_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(count_cores)
clusterExport(cl, list("prot", "params", "splitted", "model", "categories", "init_method"), envir = environment())
clusterCall(cl, function() library(evolike))

if (diff_cats){
  em_results_list <- parLapply(cl, cats, function(category){
    sink (file.path(getwd(), "output","wood_likelihood", model, paste(c(prot, "_", init_method, "_", category), collapse=""),fsep = .Platform$file.sep))
    em_results <- em_procedure(data=splitted, params=params, model = model, iter = 1000, cluster.number= category, init_method = init_method, mutation_position = "middle",  filtering = "single", trace = FALSE)
    sink() 
    em_results
  })
} else {
  em_results_list <- parLapply(cl, seq(1, trials, 1), function(trial){ 
    sink (file.path(getwd(), "output","wood_likelihood", model, paste(c(prot, "_", init_method, "_", categories, "_", trial), collapse=""),fsep = .Platform$file.sep))
    em_results <- em_procedure(data=splitted, params=params, model = model, iter = 1000, cluster.number= categories, init_method = init_method, mutation_position = "middle",  filtering = "single", trace = FALSE)
    sink() 
    em_results
  })
}

stopCluster(cl)

 