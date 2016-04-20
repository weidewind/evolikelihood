## EM algorithm


## prepare data


## clusterize observed MLE parameters
## construct rk vectors (1/0) for EM (based on cluster membership)
## find initial a (and b) parameters for EM (maximisation)
## compute initial weights based on (1/0) rk vectors (expectation)
## returns list with two values: vector of weights and matrix of parameters

params_list_to_df <- function(params){
  if (class(params) != "list"){
    stop (paste (c("Invalid argument type: ", class(params), "instead of list" )))
  }
  params <- data.frame(matrix(unlist(params), nrow=length(params), byrow=T),stringsAsFactors=FALSE)
  names(params) <- c("node", "lambda_exp_root", "lambda_root", "p_root", "p_precision" )
  params <- transform(params, lambda_exp_root = as.numeric(lambda_exp_root), lambda_root = as.numeric(lambda_root), p_root = as.numeric(p_root), p_precision = as.numeric(p_precision))
  
}


initialize_by_clustering <- function (data, parameters, mutation_position = "middle", cluster.number = 4){
  params <- params_list_to_df(parameters)
  filtered <-params[!is.na(params$p_precision) ,]
  filtered <-filtered[filtered$p_precision< 1e-5  ,]
  df <-filtered[,c("lambda_root", "p_root")]
  
  clusters <-kmeans(df, cluster.number)
  categories <- seq(from = 1, to = cluster.number, by = 1)
  rkvectors <- sapply(categories, function(k){
    rk <-sapply(clusters$cluster, function(category){
      if (category == k){1} 
      else {0}
    })
  })
  rownames(rkvectors) <- filtered$node
  
  parameters(data, mutation_position = "end", rkvector, filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
  

  
  
  
  list(iparameters = parameters, iweights = weights)
}

## alternative: choose random a (and b) 
## set equal weights for all categories
## returns list with two values: vector of weights and matrix of parameters

initialize_random <- function (){

  list(parameters = parameters, weights = weights)
}


## initialization interface
initialize <- function (method = "cluster"){
  if (method == "cluster") {
    initialize_by_clustering()
  }
  else if (method == "random"){
    initialize_random()
  }
  else {
    stop (paste ("Invalid initialization method ", method ))
  }
  
}

## EM: E - compute rk vectors and weights of each category
##     M - given rk, compute new sets of parameters for each category










### Procedures 
prot <- "h1"
prot_data <-  read.csv(paste(c(getwd(), "/input/" ,prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
params <-parameters(splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
params <- data.frame(matrix(unlist(params), nrow=length(params), byrow=T),stringsAsFactors=FALSE)
names(params) <- c("node", "lambda_exp_root", "lambda_root", "p_root", "p_precision" )
params <- transform(params, lambda_exp_root = as.numeric(lambda_exp_root), lambda_root = as.numeric(lambda_root), p_root = as.numeric(p_root), p_precision = as.numeric(p_precision))
filtered <-params[!is.na(params$p_precision) ,]
filtered <-filtered[filtered$p_precision< 1e-5  ,]
#filtered <-filtered[filtered$p_root< 30  ,]
df <-filtered[,c("lambda_root", "p_root")]
plot(df$p_root, df$lambda_root, main = "n2")
plot(df$p_root, df$lambda_root, xlim = c(0, 1.5), ylim = c(0, 0.1), main = "n2")
h1_kmeans <-kmeans(df, 3)
