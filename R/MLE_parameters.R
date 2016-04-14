install.packages("rootSolve")
library("rootSolve")
install.packages("rbenchmark")
library("rbenchmark")


## multiroot, finds only one root 
find_single_root <- function(data, mutation_position, rkvector, jack = FALSE, verbose=TRUE){
  # data is a list of dataframes, node_data is a dataframe

  if (length(data) == 1){
    mode <- "single"
    node_data <- data[[1]]
    rkvector <- c(1)
    names(rkvector) <- names(data)
    if (verbose){
      print ("Running in single mode")
    }
  }
  else {
    mode <- "group"
    if (verbose){
      print ("Running in group mode")
    }
    if(is.null(rkvector)){
      stop ("In group mode argument rkvector must be provided")
    }
  }
  
  
  if (mode == "group" || (mode == "single" && sum(node_data$event) != 0)){
    if (verbose && mode == "single"){
      print (node_data)
      print (paste("site ", node_data[2,1], " node ", node_data[2,2]))
      print (paste("number of mutations ", sum(node_data$event)))
    }
    
    pars <- list(data = data, rkvector, mutation_position = mutation_position)
    lambda_exp_root <- lambda_derivative_exp(pars)
    if (verbose){print(paste("expon_lambda root ", lambda_exp_root))}
    
    # if p root is negative, search for positive root
    p = -1
    init = 1
    count = 1
    
    while(!is.na(p) && p < 0 && count < 6){
      if (jack){
        solution_p <- multiroot(f = p_derivative, start = c(init), jacfunc = p_derivative_jacfunc, jactype = "fullusr", parms = pars)
      }
      else {
        solution_p <- multiroot(f = p_derivative, start = c(init), parms = pars)
      }
      p <-solution_p$root
      init <- init + 3
      count <- count + 1
    }
    
    if (verbose){print(paste("precision_p ", solution_p$estim.precis) )}
    if(is.na(solution_p$root)){
      print(paste("no p_roots found"))
      c(p_root = NA, p_precision = NA,
        lambda_root = NA, lambda_exp_root = lambda_exp_root)
    }
    else if(is.na(solution_p$estim.precis)){
      if (verbose){
        print(paste("p_root is ",solution_p$root, " but estimated precision for p is Na, won't try to estimate lambda"))
      }
      c(p_root = solution_p$root, p_precision = NA,
        lambda_root = NA, lambda_exp_root = lambda_exp_root)
    }
    
    else {
      pars <- list(p = solution_p$root, data = data, rkvector, mutation_position = mutation_position)
      lambda_root <- lambda_derivative_weib(pars)
      if (verbose){ print (c(p_root = solution_p$root, lambda_root = lambda_root))}
      c(p_root = solution_p$root, p_precision = solution_p$estim.precis,
        lambda_root = lambda_root,
        lambda_exp_root = lambda_exp_root)
    }
    
  }
  else {
    if (verbose && mode == "single"){ 
      print (paste("site ", node_data[2,1], " node ", node_data[2,2]))
      print ("No mutations in the subtree, all roots NA")
    }
    c(p_root = NA, p_precision = NA,
      lambda_root = NA, lambda_exp_root = NA)
  }
}

## computes parameters for all nodes, outputs only complete sets of parameters
#parameters <-function(prot, tag, fishy = FALSE){
parameters <-function(data, mutation_position = "end", fishy = FALSE, filter = TRUE, jack = FALSE, verbose = FALSE){
  
  ps <- lapply (names(data), function(elm, mut_pos){
    mutation_position <- mut_pos
    node_data <- data[elm]
    node_roots <- find_single_root(node_data, mutation_position, jack = jack, verbose = verbose)
    if (!filter || filter && !is.na(node_roots) &&  all(is.finite(node_roots)) && node_roots["p_precision"] < 1e-5 ) {
        c(node = elm, lambda_exp = node_roots["lambda_exp_root"], lambda_weib = node_roots["lambda_root"], p = node_roots["p_root"], p_precision = node_roots["p_precision"])
    }

  }, mut_pos = mutation_position)
  
  if (filter) {
    ps  <- Filter(Negate(is.null), ps)
  }
  ps
}

prot <- "h1"
prot_data <-  read.csv(paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/",prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)


h1_prms2 <-parameters(splitted, fishy = TRUE, filter= FALSE)
h1_prms_jack <-parameters(splitted, fishy = TRUE, jack = TRUE, filter= FALSE)
h1_prms_no_negative_roots <-parameters(splitted, fishy = TRUE, filter= FALSE)

sink("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/wtf_with_precision_h1")
h1_wtf_precision <-parameters(splitted, fishy = TRUE, filter= FALSE, verbose = TRUE)
sink()

snbenchmark(parameters(splitted, fishy = TRUE, jack = TRUE, filter= FALSE), parameters(splitted, fishy = TRUE, filter= FALSE),  replications = 1)


# how much does it take to compute all jackobians?
test_jack <- function(x, splitted){
  ps <- lapply (names(splitted), function(elm){
    node_data <- splitted[[elm]]
    parms <- list(node_data = node_data, mutation_position = "end")
    node_jacks <- p_derivative_jacfunc(x, parms)
  })
}

benchmark(parameters(splitted, fishy = TRUE, jack = TRUE, filter= FALSE), parameters(splitted, fishy = TRUE, filter= FALSE),  replications = 1)

benchmark(test_jack(1, splitted), test_jack(2, splitted), replications = 3)




lrts <-lrt_procedure(splitted, "testtag", fishy = TRUE, threshold = 3.84, parameters = h1_prms2)

##
try_df <- data.frame(matrix(unlist(h1_prms), nrow=length(h1_prms), byrow=T),stringsAsFactors=FALSE)
names(try_df) <- c("node", "lambda_exp_root", "lambda_root", "p_root", "p_precision" )
## equivalent to
h1nodes <- sapply(h1_prms, function(elm){
    elm["node"]
})
h1lambdaexp <- sapply(h1_prms, function(elm){
  as.numeric(elm["lambda_exp.lambda_exp_root"])
})
h1lambdawei <- sapply(h1_prms, function(elm){
  as.numeric(elm["lambda_weib.lambda_root"])
})
h1p <- sapply(h1_prms, function(elm){
  as.numeric(elm["p.p_root"])
})
h1pprec <-sapply(h1_prms, function(elm){
  as.numeric(elm["p_precision.p_precision"])
})

dframe <- data.frame(node = h1nodes, lambda_exp_root = h1lambdaexp,lambda_root = h1lambdawei, p_root =  h1p, p_precision = h1pprec)
##
head(dframe, 50)



## Inconvenient. computes all parameters, plots p and log(p) histograms
parameter_hists <- function(data, mutation_position = "end", fishy = FALSE){
  params <- parameters(data, mutation_position = "end", fishy)
  p_hists(params)
}



## plots p and log(p) histograms given a params list (computed by parameters function)

p_hists <- function(params){
  params <- params[sapply(params, function(elm){!is.na(elm[["p.p_root"]])})] #select nodes with p root defined
  p_roots <- sapply(params, function(elm) {
    as.numeric(elm[["p.p_root"]])
  })
  h <- hist(p_roots, breaks = 30, plot = FALSE)
  plot(h,  main = paste("Histogram of ",prot, " p"), xlab = "p")
  lh <- hist(log(p_roots), breaks = 30, plot = FALSE)
  plot(lh,  main = paste("Histogram of ",prot, " log(p)"), xlab = "log(p)")
}


## selects lamdas corresponding to  p>threshold if right == TRUE 
## ( p <= threshold if right == FALSE)
## and plots a histogram 
lambda_hist <- function(prot, params, threshold, right = TRUE){
  params <- params[sapply(params, function(elm){!is.na(elm[["p.p_root"]])})] #select nodes with p root defined
  p_roots <- sapply(params, function(elm) {
    as.numeric(elm[["p.p_root"]])
  })
  boolean <- sapply(params, function(elm) {
    
    if (as.numeric(elm[["p.p_root"]]) > threshold) {right == TRUE}
    else {right == FALSE}
  })
  
  lambda <- sapply(params[boolean], function(elm) {
    as.numeric(elm["lambda_weib.lambda_root"])
  })
  if (right){
    sign = ">"
  }
  else {
    sign = "<="
  }
  h <- hist(lambda, breaks = 30, plot = FALSE)
  plot(h,  main = paste("Histogram of ",prot, " lambda for p", sign, " ", threshold), xlab = "lambda")
  h <- hist(log(lambda), breaks = 30, plot = FALSE)
  plot(h,  main = paste("Histogram of ",prot, " log(lambda) for p", sign, " ", threshold), xlab = "lambda")
}
