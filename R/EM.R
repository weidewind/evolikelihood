## EM algorithm


## prepare data


## clusterize observed MLE parameters
## construct rk vectors (1/0) for EM (based on cluster membership)
## find initial a (and b) parameters for EM (maximisation)
## compute initial weights based on (1/0) rk vectors (expectation)
## compute initial rk vectors (expectation)
initialize_by_clustering <- function (){
  
}

## alternative: choose random a (and b), 
## set equal weights for all categories,
## compute initial rk vectors (expectation)
initialize_random <- function (){
  
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
    stop (paste ("Unknown initialization method ", method ))
  }
  
}

## EM: E - compute weights of each category and rk vectors
##     M - given rk, compute new sets of parameters for each category










### Procedures 
parameters <- n2_prms 
parameters <- data.frame(matrix(unlist(parameters), nrow=length(parameters), byrow=T),stringsAsFactors=FALSE)
names(parameters) <- c("node", "lambda_exp_root", "lambda_root", "p_root", "p_precision" )
parameters <- transform(parameters, lambda_exp_root = as.numeric(lambda_exp_root), lambda_root = as.numeric(lambda_root), p_root = as.numeric(p_root), p_precision = as.numeric(p_precision))
filtered <-parameters[!is.na(parameters$p_precision) ,]
filtered <-filtered[filtered$p_precision< 1e-5  ,]
#filtered <-filtered[filtered$p_root< 30  ,]
df <-filtered[,c("lambda_root", "p_root")]
plot(df$p_root, df$lambda_root, main = "n2")
plot(df$p_root, df$lambda_root, xlim = c(0, 1.5), ylim = c(0, 0.1), main = "n2")
h1_kmeans <-kmeans(df, 3)