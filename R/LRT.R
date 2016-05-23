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
                             talpha_middle <- (talpha1+talpha0)/2
                             log(lambda)+log(talpha1-talpha_middle) #changed
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
  #added
  if(fishy){
      apply_res3 <-  apply( alpha_branches, 1,
                        function(elm){ 
                          if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                            talpha1 <- as.numeric(elm["t_branch_end"])
                            talpha0 <- as.numeric(elm["t_branch_start"])
                            talpha_middle <- (talpha1+talpha0)/2
                            talpha0 -  talpha_middle
                          }  
                          else { 0 }
                        })
  }
  
  if (fishy){
      lnL <-sum(apply_res1) + lambda*(sum(apply_res2)+sum(apply_res3))
  }
  else{
      lnL <-sum(apply_res1) + lambda*sum(apply_res2)
  }
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

                             #log(p) + p*log(lambda)+(p-1)*log( talpha_middle)+log(talpha1-talpha0)
                             log(p) + p*log(lambda)+(p-1)*log( talpha_middle)+log(talpha1-talpha_middle)
                             #print ("p to die")
                             #print (log(p) + p*log(lambda)+(p-1)*log( talpha_middle)+log(talpha1-talpha_middle))
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
                            if (tbeta1 == 0){0}
                            else {
                              (tbeta1^p)*( ((tbeta0/tbeta1)^p) - 1 )  
                            }
                            #print ("p to survive on beta")
                            #print ((tbeta1^p)*( ((tbeta0/tbeta1)^p) - 1 ) )
                          }  
                          else { 0 }
                        })
  ## added
  if (fishy){
      apply_res3 <-  apply( alpha_branches, 1,
                        function(elm){ 
                          if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                            talpha1 <- as.numeric(elm["t_branch_end"])
                            talpha0 <- as.numeric(elm["t_branch_start"])
                            talpha_middle <- (talpha1+talpha0)/2
                            #print ("p to survive on alpha")
                            #print ((talpha_middle^p)*( ((talpha0/talpha_middle)^p) - 1 ))
                            (talpha_middle^p)*( ((talpha0/talpha_middle)^p) - 1 )  
                          }  
                          else { 0 }
                        })
  }
  
  
  if (fishy){
      lnL <-sum(apply_res1) + (lambda^p)*(sum(apply_res2)+sum(apply_res3))
  }
  else {
      lnL <-sum(apply_res1) + (lambda^p)*sum(apply_res2)
  }

  c(lnL = lnL, AIC = aic(lnL, 2))
}

# lnl - loglikelihood, p - number of parameters in the model
aic <- function (lnL, p){
  2*(p-lnL)
}

# lnl - loglikelihood, p - number of free parameters in the model, n - sample size
bic <- function (lnL, model = NULL, cluster.number, n){
  if (model == "weibull"){
    p <- 3*cluster.number -1;
  }
  else {
    p <- 2*cluster.number -1;
  }
  -2*lnL + p*log(n)
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
  
  out <- c(exponential[1], weibull[1], LR)
  names(out) = c("exponential", "weibull", "lr")
  out
}




## Compute LR
## takes an already computed set of parameters (with parameters(..., filter=FALSE))
## OR, if parameters = NULL, computes them 
## ! mutation_position and parameters cannot be defined simultaneously
## threshold - significance threshold (default - 3.84, value of chi-square distr with 1 degree of freedom for sign level 0.05)
## Prints parameters and lr for significant nodes
## outputs parameters and lr for all nodes
## files like "h1_single_root_31_03_truelambda_fishy" were printed here

lrt_procedure <-function(data, prot, tag, fishy = FALSE, threshold = 3.84,  mutation_position = "end", pack = "rootsolve", params = NULL, all = FALSE){
 # if (!is.null(params)){
 #   warning("Parameters are defined, therefore mutation_position will be ignored")
  #}
  sink(file = paste(c(getwd(), "/output/" ,prot,"_single_root_", tag), collapse=""), append = FALSE, type = c("output", "message"),
       split = FALSE)
  
  
  lratios <- lapply (names(data), function(elm, mutation_position){
    mutation_position <- mutation_position
    node_data <- data[elm]

    if (is.null(params)){
      node_roots <- as.list(find_single_root(node_data, mutation_position, pack=pack))
    }
    else {
      node_roots <- params[params$node == elm,]
    }
    if(!is.na(node_roots) && all(!is.na(node_roots)) && node_roots$p_precision < 1e-5  ){
      
      lr <-lrt (node_data, lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, fishy= TRUE, verbose=FALSE)
      row <- c(node = elm,  lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, p_precision = node_roots$p_precision, exponential =lr["exponential"], weibull=lr["weibull"], lr = lr["lr"])
      if(!is.na(lr["lr"]) && lr["lr"] > threshold){
        if (all){
          print (c(node = elm,  lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, p_precision = node_roots$p_precision, exponential =lr["exponential"], weibull=lr["weibull"], lr = lr["lr"]))
          
        } else print (c(node = elm,  lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, p_precision = node_roots$p_precision, lr = lr["lr"]))
      }
    } 
    else {
      row <- c(node = elm,  lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, p_precision = node_roots$p_precision, lr = NA) 
    }
  }, mutation_position)
  sink()
  
  lratios
}


dumb_wood_likelihood <-function(data, prot, tag, fishy = FALSE, threshold = 0,  mutation_position = "end", pack = "rootsolve", params, all = FALSE){

  sink(file = paste(c(getwd(), "/output/" ,prot,"_dumb_wood_", tag), collapse=""), append = FALSE, type = c("output", "message"),
       split = FALSE)
  
  weib_sum <- 0
  exp_sum <- 0
  for (elm in names(data)){
    node_data <- data[elm]
    node_roots <- params[params$node == elm,]
    
    if(nrow(node_roots) > 0 && !is.na(node_roots) && all(!is.na(node_roots)) && node_roots$p_precision < 1e-5  ){
      lr <-lrt (node_data, lambda_exp = node_roots$lambda_exp_root, lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, fishy= TRUE, verbose=FALSE)
      print(lr)
      weib_sum <- weib_sum + lr["weibull"]
      exp_sum <- exp_sum + lr["exponential"]
    }
  }
  print(exp_sum) 
  print (weib_sum)
  sink()

}




## verbose mode also prints a file with parameters for all nodes

lrt_all <- function(mutation_position = "middle", fishy = TRUE,  pack = "rootsolve", tag = "defaulttag", verbose = TRUE, threshold = 3.84, all = FALSE){
  prots <- c("h1", "h3", "n1", "n2")
  for (pr in prots){
    prot <- pr
    prot_data <-  read.csv(paste(c(getwd(), "/data/" ,prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
    splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
    if (verbose){
      sink(file = paste(c(getwd(), "/output/" ,prot,"_all_parms_", tag), collapse=""))
    }
    prms <-parameters(splitted,  mutation_position = mutation_position,  verbose = verbose, pack = pack, filter= FALSE)
    if (verbose){sink()}
    lrt_procedure(data = splitted, prot = prot, tag = tag, fishy=fishy, params = prms, threshold = threshold, all = all)
    sink()
  }
}

###
## tests and procedures


#lrt_all(mutation_position = "middle", fishy = TRUE, tag = "middle_search", pack = "rootsolve", verbose = TRUE)


#benchmark(parameters(splitted,  jack = FALSE, pack = "nleqslv", filter= FALSE), parameters(splitted,  jack = FALSE, pack = "rootsolve", filter= FALSE),  replications = 1)

##check that there are no zero-length branches with mutations
no_muts_on_zero_branches <-function(prot){
  prot_data <-  read.csv(paste(c(getwd(), "/data/" ,prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
  nullength <- prot_data[prot_data["event_indicator"] == 1,]
  nullength <-nullength[nullength["t_branch_end"]-nullength["t_branch_start"] == 0,]
  if (nrow(nullength) == 0){TRUE}
  else {FALSE}
}

# how much does it take to compute all jackobians?
test_jack <- function(x, splitted){
  ps <- lapply (names(splitted), function(elm){
    node_data <- splitted[[elm]]
    parms <- list(node_data = node_data, mutation_position = "end")
    node_jacks <- p_derivative_jacfunc(x, parms)
  })
}

#prot <- "h1"
#prot_data <-  read.csv(paste(c(getwd(), "/input/" ,prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
#splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
#node_data <- splitted["151.INTNODE4195"]
#node_roots <- find_single_root(node_data, mutation_position = "middle", jack = FALSE, pack = "nleqslv", verbose = TRUE)

#draw_hazard <-function(data=splitted, nodename = "78.NTNODE4232", to = 20, by = 0.01, mutation_position = "middle", fishy = TRUE)
#draw_lnlikelihood (data=splitted, nodename = "169.INTNODE2065", to = 20, by = 0.01, mutation_position = "middle", fishy = TRUE)
#draw_hazard(data=splitted, nodename = "78.INTNODE4232", to = 20, by = 0.05, mutation_position = "middle", fishy = TRUE)
#h1_params <-parameters(splitted, mutation_position = "middle", filter= FALSE)
#h1_prms_jack <-parameters(splitted, fishy = TRUE, jack = TRUE, filter= FALSE)
#h1_prms_no_negative_roots <-parameters(splitted, fishy = TRUE, filter= FALSE)

#sink("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/likelihood_games_h1")
#h1_likelihood_games <-parameters(splitted, mutation_position = "middle", filter= FALSE, verbose = TRUE)
#sink()
#lrt_procedure(data = splitted, prot = prot, tag = "likelihood_games", fishy=TRUE, params = h1_likelihood_games)
#sink()
