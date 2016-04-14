## fishy = FALSE: ln of true likelihood (mutation happens somewhere at the branch, no presumptions about where exactly)
## of the node_data given exponential model with parameter lambda
## fishy = TRUE: ln of likelihood computed according to eq (1) on page 1,
## Adot taken at the middle of the branch, (talpha1-talpha0) = length of the branch
lnlikelihood_exp <-function(node_data, lambda, fishy = FALSE){
  if (class(node_data) == "list"){
    node_data <- node_data[[1]]
  }
  beta_branches = node_data[node_data["event_indicator"]==0,]
  alpha_branches = node_data[node_data["event_indicator"]==1,]
  if (fishy){
    apply_res1 <- apply( alpha_branches, 1,
                         function(elm){ 
                           if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                             talpha1 <- as.numeric(elm["t_branch_end"])
                             talpha0 <- as.numeric(elm["t_branch_start"])
                             log(lambda)+log(talpha1-talpha0)
                           }  
                           else { 0 }
                         })
  }
  else {
    apply_res1 <- apply( alpha_branches, 1,
                         function(elm){ 
                           if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                             talpha1 <- as.numeric(elm["t_branch_end"])
                             talpha0 <- as.numeric(elm["t_branch_start"])
                             survival <- exp(lambda*(talpha0-talpha1))
                             log(1-survival)
                           }  
                           else { 0 }
                         })
  }
  apply_res2 <-  apply( beta_branches, 1,
                        function(elm){ 
                          if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                            tbeta1 <- as.numeric(elm["t_branch_end"])
                            tbeta0 <- as.numeric(elm["t_branch_start"])
                            tbeta0 -  tbeta1
                          }  
                          else { 0 }
                        })
  
  lnL <-sum(apply_res1) + lambda*sum(apply_res2)
  c(lnL = lnL, AIC = aic(lnL, 1))
}

## fishy = FALSE: ln of true likelihood (mutation happens somewhere at the branch, no presumptions about where exactly)
## of the node_data given weibull model with parameters lambda and p
## fishy = TRUE: ln of likelihood computed according to eq (1) on page 1,
## Adot taken at the middle of the branch, (talpha1-talpha0) = length of the branch
lnlikelihood_weibull <-function(node_data, lambda, p, fishy = FALSE){
  if (class(node_data) == "list"){
    node_data <- node_data[[1]]
  }
  beta_branches = node_data[node_data["event_indicator"]==0,]
  alpha_branches = node_data[node_data["event_indicator"]==1,]
  if (fishy){
    apply_res1 <- apply( alpha_branches, 1,
                         function(elm){ 
                           if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                             talpha1 <- as.numeric(elm["t_branch_end"])
                             talpha0 <- as.numeric(elm["t_branch_start"])
                             talpha_middle <- (talpha1+talpha0)/2
                             log(p) + p*log(lambda)+(p-1)*log(talpha_middle)+log(talpha1-talpha0)
                           }  
                           else { 0 }
                         })
    
  }
  else {
    apply_res1 <- apply( alpha_branches, 1,
                         function(elm){ 
                           if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                             talpha1 <- as.numeric(elm["t_branch_end"])
                             talpha0 <- as.numeric(elm["t_branch_start"])
                             survival <- exp((lambda^p)*(talpha1^p)*( ((talpha0/talpha1)^p) - 1 ))
                             log(1-survival)
                           }  
                           else { 0 }
                         })
  }
  apply_res2 <-  apply( beta_branches, 1,
                        function(elm){ 
                          if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                            tbeta1 <- as.numeric(elm["t_branch_end"])
                            tbeta0 <- as.numeric(elm["t_branch_start"])
                            (tbeta1^p)*( ((tbeta0/tbeta1)^p) - 1 )  
                          }  
                          else { 0 }
                        })
  
  lnL <-sum(apply_res1) + (lambda^p)*sum(apply_res2)
  c(lnL = lnL, AIC = aic(lnL, 2))
}

# lnl - loglikelihood, p - number of parameters in the model
aic <- function (lnL, p){
  2*(p-lnL)
}

# LR  LR = 2*(lnL1-lnL2) approximately follows a chi-square distribution with 1 degree of freedom
lrt <- function (node_data, lambda_exp, lambda_weib, p, fishy = FALSE, verbose = FALSE){
  weibull <- lnlikelihood_weibull(node_data, lambda = lambda_weib, p = p, fishy)
  exponential <- lnlikelihood_exp(node_data, lambda = lambda_exp, fishy)
  LR <- 2*(weibull[1]-exponential[1])
  if(verbose){
    print ("weibull")
    print (weibull)
    print ("exponential")
    print (exponential)
    print ("LR")
    print (LR)
  }
  LR
}




## Compute LR
## takes an already computed set of parameters (with parameters(..., filter=FALSE))
## OR, if parameters = NULL, computes them 
## ! mutation_position and parameters cannot be defined simultaneously
## threshold - significance threshold (default - 3.84, value of chi-square distr with 1 degree of freedom for sign level 0.05)
## Prints parameters and lr for significant nodes
## outputs parameters and lr for all nodes
## files like "h1_single_root_31_03_truelambda_fishy" were printed here

lrt_procedure <-function(data, prot, tag, fishy = FALSE, threshold = 3.84,  mutation_position = "end", parameters = NULL){
  if (!is.null(parameters)){
    warning("Parameters are defined, therefore mutation_position will be ignored")
  }
  sink(file = paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/",prot,"_single_root_", tag), collapse=""), append = FALSE, type = c("output", "message"),
       split = FALSE)
  
  if (!is.null(parameters)){
    parameters <- data.frame(matrix(unlist(parameters), nrow=length(parameters), byrow=T),stringsAsFactors=FALSE)
    names(parameters) <- c("node", "lambda_exp_root", "lambda_root", "p_root", "p_precision" )
    parameters <- transform(parameters, lambda_exp_root = as.numeric(lambda_exp_root), lambda_root = as.numeric(lambda_root), p_root = as.numeric(p_root), p_precision = as.numeric(p_precision))
  }
  
  lratios <- lapply (names(data), function(elm, mutation_position){
    mutation_position <- mutation_position
    node_data <- data[elm]
    if (is.null(parameters)){
      node_roots <- as.list(find_single_root(node_data, mutation_position))
    }
    else {
      node_roots <- parameters[parameters$node == elm,]
    }
    if(!is.na(node_roots) && all(!is.na(node_roots)) && node_roots$p_precision < 1e-5  ){
      
      lr <-lrt (node_data, lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_root, p = node_roots$p_root, fishy= TRUE, verbose=FALSE)
      row <- c(node = elm,  lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_root, p = node_roots$p_root, p_precision = node_roots$p_precision, lr = lr)
      if(!is.na(lr) && lr > threshold){
        print (c(node = elm,  lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_root, p = node_roots$p_root, p_precision = node_roots$p_precision, lr = lr))
      }
    } 
    else {
      row <- c(node = elm,  lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_root, p = node_roots$p_root, p_precision = node_roots$p_precision, lr = NA) 
    }
  }, mutation_position)
  sink()
  
  lratios
}

prot <- "h1"
prot_data <-  read.csv(paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/",prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
h1_prms <-parameters(splitted, fishy = TRUE, filter= FALSE)
lrt_procedure(data = splitted, prot = prot, tag = "without_negative_roots_step_3", fishy=TRUE, parameters = h1_prms)
sink()
prot <- "h3"
prot_data <-  read.csv(paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/",prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
h3_prms <-parameters(splitted, fishy = TRUE, filter= FALSE)
lrt_procedure(data = splitted, prot = prot, tag = "without_negative_roots", fishy=TRUE, parameters = h3_prms)
sink()
prot <- "n1"
prot_data <-  read.csv(paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/",prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
n1_prms <-parameters(splitted, fishy = TRUE, filter= FALSE)
lrt_procedure(data = splitted, prot = prot, tag = "without_negative_roots", fishy=TRUE, parameters = n1_prms)
sink()
prot <- "n2"
prot_data <-  read.csv(paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/",prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
n2_prms <-parameters(splitted, fishy = TRUE, filter= FALSE)
lrt_procedure(data = splitted, prot = prot, tag = "without_negative_roots", fishy=TRUE, parameters = n2_prms)
sink()


###


