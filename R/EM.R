list.of.packages <- c("rgl")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages)

## EM algorithm


## prepare data


## clusterize observed MLE parameters
## construct rk vectors (1/0) for EM (based on cluster membership)
## find initial a (and b) parameters for EM (maximisation)
## compute initial weights based on (1/0) rk vectors (expectation)
## returns list with two values: vector of weights and matrix of parameters



initialize_by_clustering <- function (data, params, mutation_position = "middle", cluster.number = 4){
#  f <- filter_unsolved(data=data, params=params)
#  fparams = f$fparams
 # fdata = f$fdata
  parclust <-params[,c("lambda_weib_root", "p_root")]
  
  #filtered <-params[!is.na(params$p_precision) ,]
 # filtered <-filtered[filtered$p_precision< 1e-5  ,]
  #df <-filtered[,c("lambda_weib_root", "p_root")]
 # data_filtered <- data[filtered$node]
  
  scaling <-max(parclust[,2])/max(parclust[,1])
  parclust_scaled <- data.frame(lambda_weib_root = scaling*parclust[,1], p_root = parclust[,2])
  
  clusters <-kmeans(parclust_scaled, cluster.number)
  categories <- seq(from = 1, to = cluster.number, by = 1)
  rkvectors <- sapply(categories, function(k){
    rk <-sapply(clusters$cluster, function(category){
      if (category == k){1} 
      else {0}
    })
  })
  rownames(rkvectors) <- params$node
  
  init_params <-compute_params(data = data, rkvectors = rkvectors, mutation_position = mutation_position)
  
  init_weights <- compute_weights(rkvectors)
  
  list(iparameters = init_params, iweights = init_weights)
}


## alternative: choose random a (and b) 
## set equal weights for all categories
## returns list with two values: vector of weights and matrix of parameters

initialize_random <- function (){
  
}

initialize_by <- function (init_params, init_weights, cluster.number = 4){
  if (class(init_params) != "numeric") {
    stop ("Invalid argument params: expected vector of length 2*cluster.number")
  }
  if (length(init_params) != 2*cluster.number) {
    stop ("Invalid params length: expected vector of length 2*cluster.number")
  }
  if (length(init_weights) != cluster.number) {
    stop ("Invalid weights length: expected vector of length cluster.number")
  }
  
  categories <- seq(from = 1, to = cluster.number, by = 1)
  init_params <-matrix(init_params, ncol =2, nrow = 4, byrow = TRUE)
  colnames(init_params) = c("lambda", "p")

  
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


em_procedure <-function(data, params, iter = 100, cluster.number= 4, init_method = c("cluster", "random", "by"), init_params = NULL, init_weights = NULL, mutation_position = "middle",  filtering = c("single", "unsolved"), trace = TRUE){
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
    init <- initialize_by_clustering(data=fdata, params=fparams, mutation_position = mutation_position, cluster.number = cluster.number)
    iparameters <- init$iparameters
    iweights <- init$iweights 
  }
  else if (init_method == "random"){
    init <- initialize_random(data=fdata, params=fparams, mutation_position = mutation_position, cluster.number = cluster.number)
    iparameters <- init$iparameters
    iweights <- init$iweights 
  }
  else if (init_method == "by"){
    if (is.null(init_params) || is.null(init_weights)){
      stop ("Chosen inititalization method requires init_params and init_weights arguments")
    }
    init <- initialize_by(init_params, init_weights, cluster.number = cluster.number)
    iparameters <- init$iparameters
    iweights <- init$iweights 
  }
  else {
    stop (paste ("Invalid initialization method ", method ))
  }
  
  em_results <- em(fdata, iparameters, iweights, iter= iter, mutation_position = mutation_position, cluster.number = cluster.number, trace = trace)

}


em <- function(data, parameters, weights, iter = 100, cluster.number= 4, mutation_position = "middle", trace = TRUE){
  #todo: data must be filtered (before initialization?), and this new set must be used in em/ params (for every node) are not to be confused with cluster paramteres.
  #  init <- initialize (data, params=parameters, mutation_position = mutation_position, cluster.number = cluster.number, method = "cluster")
  #  parameters <- init$iparameters
  #  weights <- init$iiweights
  
  #  dev.new()
  #  h1=dev.cur()
  #  myplot1 <- scatterplot3d(parameters[1,"lambda"], parameters[1,"p"], weights[1], color= "red", type="h", pch=19, xlim = c(max(parameters[1,"lambda"] - 0.1,0),parameters[1,"lambda"] + 0.1), ylim= c(max(parameters[1,"p"] - 1,0),parameters[1,"p"] + 1), zlim = c(0,1))
  #  dev.new()
  #  h2=dev.cur()
  #  myplot2 <- scatterplot3d(parameters[2,"lambda"], parameters[2,"p"], weights[2], color= "blue", type="h", pch=19, xlim = c(max(parameters[2,"lambda"] - 0.1,0),parameters[2,"lambda"] + 0.1), ylim= c(max(parameters[2,"p"] - 1,0),parameters[2,"p"] + 1), zlim = c(0,1))
  #  dev.new()
  #  h3=dev.cur()
  # myplot3 <- scatterplot3d(parameters[3,"lambda"], parameters[3,"p"], weights[3], color= "green", type="h", pch=19, xlim = c(max(parameters[3,"lambda"] - 0.1,0),parameters[3,"lambda"] + 0.1), ylim= c(max(parameters[3,"p"] - 1,0),parameters[3,"p"] + 1), zlim = c(0,1))
  # dev.new()
  #  h4=dev.cur()
  # myplot4 <- scatterplot3d(parameters[3,"lambda"], parameters[4,"p"], weights[4], color= "black", type="h", pch=19, xlim = c(max(parameters[4,"lambda"] - 0.1,0),parameters[4,"lambda"] + 0.1), ylim= c(max(parameters[4,"p"] - 1,0),parameters[4,"p"] + 1), zlim = c(0,1))
  # dev.new()
  # h5=dev.cur()
  if (trace){
   myplot <- scatterplot3d(parameters[1,"lambda"], parameters[1,"p"], weights[1], color= "red", type="h", xlim = c(0, 0.1), ylim = c(0,10), zlim = c(0,1), pch=19)
    # dev.set(h5)
    myplot$points3d(parameters[2,"lambda"], parameters[2,"p"], weights[2], col= "blue", pch=19, type="h")
    myplot$points3d(parameters[3,"lambda"], parameters[3,"p"], weights[3], col= "green", pch=19, type="h")
    myplot$points3d(parameters[4,"lambda"], parameters[4,"p"], weights[4], col= "black",pch=19,  type="h")
  }
  old_lnL <- NULL
  for (i in seq(1,iter,1)){
    print (paste(c("------------Step ", i), collapse=""))
    rkvectors <- compute_rkvectors(data=data, parameters=parameters, weights=weights)
    parameters <- compute_params(data=data, rkvectors=rkvectors, mutation_position = mutation_position)
    weights <- compute_weights(rkvectors)
    print(rkvectors)
    print(weights)
    # dev.set(h1)
    #  myplot1$points3d(parameters[1,"lambda"], parameters[1,"p"], weights[1], col= "red", type="h")
    # dev.set(h2)
    # myplot2$points3d(parameters[2,"lambda"], parameters[2,"p"], weights[2], col= "blue", type="h")
    # dev.set(h3)
    #  myplot3$points3d(parameters[3,"lambda"], parameters[3,"p"], weights[3], col= "green", type="h")
    # dev.set(h4)
    # myplot4$points3d(parameters[4,"lambda"], parameters[4,"p"], weights[4], col= "black", type="h")
    
    # dev.set(h5)
    if (trace){
      myplot$points3d(parameters[1,"lambda"], parameters[1,"p"], weights[1], col= "red", type="h")
      myplot$points3d(parameters[2,"lambda"], parameters[2,"p"], weights[2], col= "blue", type="h")
      myplot$points3d(parameters[3,"lambda"], parameters[3,"p"], weights[3], col= "green", type="h")
      myplot$points3d(parameters[4,"lambda"], parameters[4,"p"], weights[4], col= "black", type="h")
    }  
    model_lnL <- compute_model_lnL(data=data, parameters=parameters, weights=weights)
    print ("model lnL")
    print(model_lnL)
    if (!is.null(old_lnL) && model_lnL - old_lnL < 0.0001){
      break
    }
    else {old_lnL <- model_lnL}
  }
  list(parameters=parameters, rkvectors=rkvectors, weights=weights, lnL = model_lnL)
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

compute_rkvectors <- function(data, parameters, weights){
  cluster.number = length(weights)
  categories <- seq(from = 1, to = cluster.number, by = 1)
  rkvectors <- sapply(categories, function(k){
    rk <-sapply(names(data), function(nodename){
      cat_probs <- sapply ( categories, function (cat) {
        lnL_dat <- lnlikelihood_weibull(data[[nodename]], parameters[cat,"lambda"], parameters[cat,"p"], fishy = TRUE)
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

compute_params <- function(data, rkvectors, mutation_position = "middle" ){
  cluster.number = ncol(rkvectors)
  categories <- seq(from = 1, to = cluster.number, by = 1)
  new_params <-matrix(nrow = cluster.number, ncol = 2)
  colnames(new_params) = c("lambda", "p")
  for (k in categories){
    k_params <- find_single_root(data = data, mutation_position=mutation_position, rkvector = rkvectors[, k], jack = FALSE, pack = "rootsolve", verbose=TRUE)
    new_params[k, "p"]  <- k_params["p_root"]
    new_params[k, "lambda"]  <- k_params["lambda_weib_root"]
  }
  new_params
}


compute_model_lnL <- function(data, parameters, weights){
  cluster.number = length(weights)
  categories <- seq(from = 1, to = cluster.number, by = 1)

    likelihood_vector <-sapply(names(data), function(nodename){
      cat_probs <- sapply ( categories, function (cat) {
        lnL_dat <- lnlikelihood_weibull(data[[nodename]], parameters[cat,"lambda"], parameters[cat,"p"], fishy = TRUE)
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
