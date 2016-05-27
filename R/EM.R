list.of.packages <- c("scatterplot3d", "parallel")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

#Don’t use library() or require(). These modify the search path, affecting what functions are available from the global environment. 
#It’s better to use the DESCRIPTION to specify your package’s requirements
#library(scatterplot3d)
## EM algorithm


## prepare data


## clusterize observed MLE parameters
## construct rk vectors (1/0) for EM (based on cluster membership)
## find initial a (and b) parameters for EM (maximisation)
## compute initial weights based on (1/0) rk vectors (expectation)
## returns list with two values: vector of weights and matrix of parameters



initialize_by_clustering <- function (data, params, mutation_position = "middle", cluster.number = 4, model = NULL){
  if (model == "weibull"){
    parclust <-params[,c("lambda_weib_root", "p_root")]
    scaling <-max(parclust[,2])/max(parclust[,1])
    parclust_scaled <- data.frame(lambda_weib_root = scaling*parclust[,1], p_root = parclust[,2])
  }
  else {
    parclust_scaled <- data.frame(lambda_exp_root = params[,c("lambda_exp_root")])
  }
  
  clusters <-kmeans(parclust_scaled, cluster.number)
  categories <- seq(from = 1, to = cluster.number, by = 1)
  rkvectors <- sapply(categories, function(k){
    rk <-sapply(clusters$cluster, function(category){
      if (category == k){1} 
      else {0}
    })
  })
  rownames(rkvectors) <- params$node
  
  init_params <-compute_params(data = data, rkvectors = rkvectors, model = model, mutation_position = mutation_position)
  
  init_weights <- compute_weights(rkvectors)
  
  list(iparameters = init_params, iweights = init_weights)
}


## alternative: choose random a (and b) 
## set equal weights for all categories
## returns list with two values: vector of weights and matrix of parameters

initialize_random <- function (params, cluster.number = 4, model = NULL){
  if (model == "weibull"){
    init_params <-matrix(ncol =2, nrow = cluster.number, byrow = TRUE)
    colnames(init_params) = c("lambda", "p")
    p_roots <- params[, "p_root"]
    prand <- exp(runif(cluster.number, min=min(log (p_roots)),max=max(log (p_roots))))
    init_params[, "p"] <-prand
    lambda_roots <- params[, "lambda_weib_root"]
    lambdarand <- exp(runif(cluster.number, min=min(log (lambda_roots)),max=max(log (lambda_roots))))
    init_params[, "lambda"] <-lambdarand
  }
  else if (model == "exponential"){
    init_params <-matrix( ncol =1, nrow = cluster.number, byrow = TRUE)
    colnames(init_params) = c("lambda")
    lambda_roots <- params[, "lambda_exp_root"]
    lambdarand <- exp(runif(cluster.number, min=min(log (lambda_roots)),max=max(log (lambda_roots))))
    init_params[, "lambda"] <-lambdarand
  }
  
  init_weights = rep(1/cluster.number, cluster.number)
  list(iparameters = init_params, iweights = init_weights)
}

initialize_by <- function (init_params, init_weights, model = NULL, cluster.number = 4){
  if (class(init_params) != "numeric") {
    stop ("Invalid argument params: expected vector of length 2*cluster.number")
  }
  if (model == "weibull"){
    if (length(init_params) != 2*cluster.number) {
      stop ("Invalid params length: expected vector of length 2*cluster.number")
    }
  }
  else if (model == "exponential"){
    if (length(init_params) != cluster.number) {
      stop ("Invalid params length: expected vector of length cluster.number")
    }
  }
  else {
    stop (paste (c("Invalid model name: expected 'weibull' or 'exponential', recieved ", model), collpase = ""))
  }
  if (length(init_weights) != cluster.number) {
    stop ("Invalid weights length: expected vector of length cluster.number")
  }
  if (sum(init_weights) != 1){
    stop ("Invalid weights value: sum of weights must be equal to 1")
  }
  
  categories <- seq(from = 1, to = cluster.number, by = 1)
  if (model == "weibull"){
    init_params <-matrix(init_params, ncol =2, nrow = cluster.number, byrow = TRUE)
    colnames(init_params) = c("lambda", "p")
  }
  else if (model == "exponential"){
    init_params <-matrix(init_params, ncol =1, nrow = cluster.number, byrow = TRUE)
    colnames(init_params) = c("lambda")
  }
  
  list(iparameters = init_params, iweights = init_weights)
}

filter_unsolved <-function(data, params){
  fparams <-params[!is.na(params$p_precision) ,]
  fparams <-fparams[fparams$p_precision< 1e-5  ,]
  fdata <- data[fparams$node]
  list(fdata = fdata, fparams = fparams)
}

filter_unsolved_and_single <-function(data, params){
  fparams <-params[!is.na(params$p_precision),]
  fparams <-fparams[fparams$p_precision< 1e-5,]
  fparams <-fparams[fparams$events > 1,]
  fdata <- data[fparams$node]
  list(fdata = fdata, fparams = fparams)
}


em_procedure <-function(data, params, model = NULL, iter = 100, cluster.number= 4, init_method = c("cluster", "random", "by"), init_params = NULL, init_weights = NULL, mutation_position = "middle",  filtering = c("single", "unsolved"), trace = TRUE){
  if (filtering == "single"){
    fi <- filter_unsolved_and_single(data=data, params=params)
    fdata <- fi$fdata
    fparams <- fi$fparams
  }
  else if (filtering == "unsolved"){
    fi <- filter_unsolved(data=data, params=params)
    fdata <- fi$fdata
    fparams <- fi$fparams
  }
  else {
    stop("Invalid argument filtering: must be either single or unsolved")
  }
  
  if (init_method == "cluster") {
    init <- initialize_by_clustering(data=fdata, params=fparams, model = model, mutation_position = mutation_position, cluster.number = cluster.number)
    iparameters <- init$iparameters
    iweights <- init$iweights 
  }
  else if (init_method == "random"){
    init <- initialize_random(params=fparams, model = model, cluster.number = cluster.number)
    iparameters <- init$iparameters
    iweights <- init$iweights 
  }
  else if (init_method == "by"){
    if (is.null(init_params) || is.null(init_weights)){
      stop ("Chosen inititalization method requires init_params and init_weights arguments")
    }
    init <- initialize_by(init_params, init_weights, model = model, cluster.number = cluster.number)
    iparameters <- init$iparameters
    iweights <- init$iweights 
  }
  else {
    stop (paste ("Invalid initialization method ", method ))
  }
  
  print ("Initial parameters:")
  print (iparameters)
  print (iweights)
  
  em_results <- em(data = fdata, parameters = iparameters, model = model, weights = iweights, iter= iter, mutation_position = mutation_position, cluster.number = cluster.number, trace = trace)

}


em <- function(data, model = NULL, parameters, weights, iter = 100, cluster.number= 4, mutation_position = "middle", trace = TRUE){
  if (trace){
    myplot <- tracer(parameters, weights, cluster.number, init = TRUE)
  }
  old_lnL <- NULL
  for (i in seq(1,iter,1)){
    print (paste(c("------------Step ", i), collapse=""))
    rkvectors <- compute_rkvectors(data=data, parameters=parameters, model = model, weights=weights)
    parameters <- compute_params(data=data, rkvectors=rkvectors, model = model, mutation_position = mutation_position)
    weights <- compute_weights(rkvectors)
    print(rkvectors)
    print(weights)
    if (trace){
      myplot <- tracer(parameters, weights, cluster.number, myplot = myplot, init = FALSE)
    }  
    model_lnL <- compute_model_lnL(data=data, model = model,parameters=parameters, weights=weights)
    print ("model lnL")
    print(model_lnL)
    if (!is.null(old_lnL) && model_lnL - old_lnL < 0.0001){
      break
    }
    else {old_lnL <- model_lnL}
  }
  
  model_bic <- bic(lnL = model_lnL, model = model, cluster.number = cluster.number, n = length(data))
  print ("data length")
  print (length(data))
  print ("model bic")
  print(model_bic)
  list(parameters=parameters, rkvectors=rkvectors, weights=weights, lnL = model_lnL, bic = model_bic)
}

tracer <- function (parameters, weights, cluster.number, myplot,  init = FALSE){
  colors <- c("red", "blue", "green", "black", "orange", "gray", "violet")
  if (ncol(parameters) == 2){ #weibull model
    if (init){
      myplot <- scatterplot3d(parameters[1,"lambda"], parameters[1,"p"], weights[1], color= colors[1], type="h", xlim = c(0, 0.1), ylim = c(0,10), zlim = c(0,1), pch=19)
      if (cluster.number > 1){
        for (i in seq(2,cluster.number,1)){
          myplot$points3d(parameters[i,"lambda"], parameters[i,"p"], weights[i], col= colors[i], pch=19, type="h")
        }
      }
    }
    else {
      for (i in seq(1,cluster.number,1)){
        print ("weights")
        print(weights)
        myplot <- myplot
        myplot$points3d(parameters[i,"lambda"], parameters[i,"p"], weights[i], col= colors[i], type="h")
      }
    }
  }
  
  else if (ncol(parameters) == 1){ #exponential model
    if (init){
      myplot <- plot(parameters[1,"lambda"], weights[1], col= colors[1],  xlim = c(0, 0.1), ylim = c(0,1), xlab = "lambda", ylab = "weight", pch=19)
      if (cluster.number > 1){
        for (i in seq(2,cluster.number,1)){
          points(parameters[i,"lambda"], weights[i], col= colors[i], pch=19)
        }
      }
    }
    else {
      for (i in seq(1,cluster.number,1)){
        points(parameters[i,"lambda"],  weights[i], col= colors[i])
      }
    }
  }
  myplot
}


compute_weights <- function(rkvectors){
  if (class(rkvectors) != "matrix"){
    stop (paste (c("Invalid type of rkvectors: expected matrix, got ", class(rkvectors))))
  }
  categories <- seq(from = 1, to = ncol(rkvectors), by = 1)
  weights <- sapply (categories, function (k){
    rkvector = rkvectors[, k]
    sum(rkvector)/length(rkvector)
  })
}

compute_rkvectors <- function(data, model = NULL, parameters, weights){
  cluster.number = length(weights)
  categories <- seq(from = 1, to = cluster.number, by = 1)
  rkvectors <- sapply(categories, function(k){
    rk <-sapply(names(data), function(nodename){
      cat_probs <- sapply ( categories, function (cat) {
        if (model == "weibull"){
          lnL_dat <- lnlikelihood_weibull(data[[nodename]], parameters[cat,"lambda"], parameters[cat,"p"], fishy = TRUE)
        }
        else {
          lnL_dat <- lnlikelihood_exp(data[[nodename]], parameters[cat,"lambda"], fishy = TRUE)
        }
        lnL <- lnL_dat[1]
        weights[cat] * exp(lnL)
      }
      )
      cat_probs[k]/sum(cat_probs)
    })
  })
  rownames(rkvectors) <- names(data)
  rkvectors
}



compute_params_insane <- function(data, model = NULL, rkvectors, mutation_position = "middle", parallel = FALSE ){
  cluster.number = ncol(rkvectors)
  categories <- seq(from = 1, to = cluster.number, by = 1)
  
  if (parallel){
    if (Sys.info()["sysname"] == "Windows"){
      count_cores <- detectCores() - 1
      cl <- makeCluster(count_cores)
      clusterExport(cl, list("data", "rkvectors", "mutation_position", "model"), envir = environment())
      clusterCall(cl, function() library(evolike))
      func <-  parLapply
    } else {
      func <- mclapply
      mc.cores <- cluster.number
    }
  } else {
    func <- mysapply
    mc.cores <- 0 # mock variable
  }
  
  params_list <- func(categories, function(k){
    k_params <- find_single_root(data = data, mutation_position=mutation_position, rkvector = rkvectors[, k], jack = FALSE, pack = "rootsolve", verbose=TRUE)
    if (model == "weibull"){
      c(k_params["lambda_weib_root"], k_params["p_root"])
    } else {
      c(k_params["lambda_exp_root"])
    }
  }, mc.cores = mc.cores)
  
  if (parallel && Sys.info()["sysname"] == "Windows"){
    stopCluster(cl)
  }
  
  if(model == "weibull"){
    new_params <- matrix(unlist(params_list), ncol = 2, byrow = TRUE)
    colnames(new_params) = c("lambda", "p")
  } else {
    new_params <- matrix(unlist(params_list), ncol = 1, byrow = TRUE)
    colnames(new_params) = c("lambda")
  }
  
  new_params
}

mysapply <- function(X, FUN, mc.cores = 1){
  sapply(X, FUN)
}



compute_params <- function(data, model = NULL, rkvectors, mutation_position = "middle"){
  cluster.number = ncol(rkvectors)
  categories <- seq(from = 1, to = cluster.number, by = 1)
  
  params_list <- sapply(categories, function(k){
    k_params <- find_single_root(data = data, mutation_position=mutation_position, rkvector = rkvectors[, k], jack = FALSE, pack = "rootsolve", verbose=TRUE)
    if (model == "weibull"){
      c(k_params["lambda_weib_root"], k_params["p_root"])
    } else {
      c(k_params["lambda_exp_root"])
    }
  })
  
  if(model == "weibull"){
    new_params <- matrix(unlist(params_list), ncol = 2, byrow = TRUE)
    colnames(new_params) = c("lambda", "p")
  } else {
    new_params <- matrix(unlist(params_list), ncol = 1, byrow = TRUE)
    colnames(new_params) = c("lambda")
  }
  
  new_params
}

compute_model_lnL <- function(data, model = NULL, parameters, weights){
  cluster.number = length(weights)
  categories <- seq(from = 1, to = cluster.number, by = 1)

    likelihood_vector <-sapply(names(data), function(nodename){
      cat_probs <- sapply ( categories, function (cat) {
        if (model == "weibull"){
          lnL_dat <- lnlikelihood_weibull(data[[nodename]], parameters[cat,"lambda"], parameters[cat,"p"], fishy = TRUE)
        }
        else {
          lnL_dat <- lnlikelihood_exp(data[[nodename]], parameters[cat,"lambda"], fishy = TRUE)
        }
        lnL <- lnL_dat[1]
        weights[cat] * exp(lnL)
      }
      )
      sum(cat_probs)
    })
  
  lnL <- sum(log(likelihood_vector))

}

## EM: E - compute rk vectors and weights of each category
##     M - given rk, compute new sets of parameters for each category










### Procedures 
#prot <- "h1"
#prot_data <-  read.csv(paste(c(getwd(), "/input/" ,prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
#splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
#params <-parameters(splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)


#params <- data.frame(matrix(unlist(params), nrow=length(params), byrow=T),stringsAsFactors=FALSE)
#names(params) <- c("node", "lambda_exp_root", "lambda_weib_root", "p_root", "p_precision" )
#params <- transform(params, lambda_exp_root = as.numeric(lambda_exp_root), lambda_weib_root = as.numeric(lambda_weib_root), p_root = as.numeric(p_root), p_precision = as.numeric(p_precision))
#filtered <-params[!is.na(params$p_precision) ,]
#filtered <-filtered[filtered$p_precision< 1e-5  ,]
##filtered <-filtered[filtered$p_root< 30  ,]
#df <-filtered[,c("lambda_weib_root", "p_root")]
#plot(df$p_root, df$lambda_weib_root, main = "n2")
#plot(df$p_root, df$lambda_weib_root, xlim = c(0, 1.5), ylim = c(0, 0.1), main = "n2")
#h1_kmeans <-kmeans(df, 3)
