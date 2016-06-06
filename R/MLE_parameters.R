list.of.packages <- c("rootSolve", "nleqslv")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')


#Don’t use library() or require(). These modify the search path, affecting what functions are available from the global environment. 
#It’s better to use the DESCRIPTION to specify your package’s requirements
#library("rootSolve") 
#library("nleqslv")



#'  ML parameters of exponential and weibull models applied to trees
#'  
#' @param data list of dataframes, produced by splitting read.csv(.._for_LRT.csv)
#' @param mutation_position end, start or middle. Defines position of mutation on a branch
#' @return A named numeric vector c(p_root, p_precision,lambda_weib_root, lambda_exp_root) 
find_single_root <- function(data, mutation_position, rkvector, jack = FALSE, pack = "rootsolve", verbose=TRUE){
  # pack nleqslv
  # data is a list of dataframes, node_data is a dataframe

  if (pack != "rootsolve" && pack != "nleqslv"){
      stop ("Incorrect argument pack: must be either 'rootsolve' or 'nleqslv'")
  }
  
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
    
    pars <- list(data = data, rkvector = rkvector, mutation_position = mutation_position)
    lambda_exp_root <- lambda_derivative_exp(pars)
    if (verbose){print(paste("expon_lambda root ", lambda_exp_root))}
    
    # if p root is negative, search for positive root
    p = -1
    init = 0.2
    count = 1
    
    while(  (!is.na(p) && (p < 0  || (!is.null(solution_p$termcd) && solution_p$termcd > 2)|| (!is.null(solution_p$estim.precis) && is.nan(solution_p$estim.precis))) ) && count < 20){
      print (paste(c("staring with p = ", init)))
      if (pack == "rootsolve"){
        if (jack){
          solution_p <- multiroot(f = p_derivative, start = c(init), jacfunc = p_derivative_jacfunc, jactype = "fullusr", parms = pars)
        }
        else {
          solution_p <- multiroot(f = p_derivative, start = c(init), parms = pars)
        }
        p <-solution_p$root
      }
      else if (pack == "nleqslv"){
        #if (jack){
        #  solution_p <- multiroot(f = p_derivative, start = c(init), jacfunc = p_derivative_jacfunc, jactype = "fullusr", parms = pars)
        #}
       # else {
          solution_p <- nleqslv(fn = p_derivative, x = c(init), parms = pars)
         # print (solution_p)
        #}
        p <-solution_p$x
      }
      
      init <- init + 1
      count <- count + 1
    }
    
    if (pack == "rootsolve"){
      precision <- solution_p$estim.precis
      p_root <- solution_p$root
    }
    else {
      precision <- solution_p$fvec
      p_root <- solution_p$x
    }
    
    
    if (verbose){print(paste("precision_p ", precision) )}
    if(is.na(p_root)){
      print(paste("no p_roots found"))
      c(p_root = NA, p_precision = NA,
        lambda_weib_root = NA, lambda_exp_root = lambda_exp_root)
    }
    else if(is.na(precision)){
      if (verbose){
        print(paste("p_root is ",p_root, " but estimated precision for p is Na, won't try to estimate lambda"))
      }
      c(p_root = p_root, p_precision = NA,
        lambda_weib_root = NA, lambda_exp_root = lambda_exp_root)
    }
    
    else {
      pars <- list(p = p_root, data = data, rkvector = rkvector, mutation_position = mutation_position)
      lambda_weib_root <- lambda_derivative_weib(pars)
      if (verbose){ print (c(p_root = p_root, lambda_weib_root = lambda_weib_root))}
      c(p_root = p_root, p_precision = precision,
        lambda_weib_root = lambda_weib_root,
        lambda_exp_root = lambda_exp_root)
    }
    
  }
  else {
    if (verbose && mode == "single"){ 
      print (paste("site ", node_data[2,1], " node ", node_data[2,2]))
      print ("No mutations in the subtree, all roots NA")
    }
    c(p_root = NA, p_precision = NA,
      lambda_weib_root = NA, lambda_exp_root = NA)
  }
}


#' @export
## computes parameters for all nodes, outputs only complete sets of parameters
#parameters <-function(prot, tag, fishy = FALSE){
parameters <-function(data, mutation_position = "middle", rkvector, filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE){
  
  ps <- lapply (names(data), function(elm, mut_pos){
    mutation_position <- mut_pos
    node_data <- data[elm]
    node_roots <- find_single_root(node_data, mutation_position, rkvector, jack = jack, pack = pack, verbose = verbose)
    if (!filter || filter && !is.na(node_roots) &&  all(is.finite(node_roots)) && node_roots["p_precision"] < 1e-5 ) {
        c(node = elm, lambda_exp = node_roots["lambda_exp_root"], lambda_weib = node_roots["lambda_weib_root"], p = node_roots["p_root"], p_precision = node_roots["p_precision"], events = sum(node_data[[1]]$event))
    }

  }, mut_pos = mutation_position)
  
  if (filter) {
    ps  <- Filter(Negate(is.null), ps)
  }
  params <- data.frame(matrix(unlist(ps), nrow=length(ps), byrow=T),stringsAsFactors=FALSE)
  names(params) <- c("node", "lambda_exp_root", "lambda_weib_root", "p_root", "p_precision", "events" )
  params <- transform(params, lambda_exp_root = as.numeric(lambda_exp_root), lambda_weib_root = as.numeric(lambda_weib_root), p_root = as.numeric(p_root), p_precision = as.numeric(p_precision), events = as.numeric(events))
  
}





## plots p and log(p) histograms given a params list (computed by parameters function)
## for dataframe params 

p_hists <- function(params){
  params <-  params[!is.na(params["p_root"]),] #select nodes with p root defined
  p_roots <- params[, "p_root"]
  h <- hist(p_roots, breaks = 30, plot = FALSE)
  plot(h,  main = paste("Histogram of ",prot, " p"), xlab = "p")
  lh <- hist(log(p_roots), breaks = 30, plot = FALSE)
  plot(lh,  main = paste("Histogram of ",prot, " log(p)"), xlab = "log(p)")
}


## selects lamdas corresponding to  p>threshold if right == TRUE 
## ( p <= threshold if right == FALSE)
## and plots a histogram 

lambda_hist <- function(prot, params, threshold, right = TRUE){
  params <-  params[!is.na(params["p_root"]),] #select nodes with p root defined
  p_roots <- params[, "p_root"]
  if (right){
    lambdas <- params[params["p_root"] > threshold, "lambda_weib_root"]
    sign = ">"
  }
  else {
    lambdas <- params[params["p_root"] <= threshold, "lambda_weib_root"]
    sign = "<="
  }

  h <- hist(lambdas, breaks = 30, plot = FALSE)
  plot(h,  main = paste("Histogram of ",prot, " lambda for p", sign, " ", threshold), xlab = "lambda")
  h <- hist(log(lambdas), breaks = 30, plot = FALSE)
  plot(h,  main = paste("Histogram of ",prot, " log(lambda) for p", sign, " ", threshold), xlab = "lambda")
}
