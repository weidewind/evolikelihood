# maximum likelihood estimation for lambda (exponential distribution)

lambda_derivative_exp <- function(parms){
  node_data <- parms$node_data
  mutation_position <- parms$mutation_position
  if (mutation_position == "start"){
    no_events_branches = node_data[node_data["event_indicator"]==0,]
  }
  else if (mutation_position == "end"){
    no_events_branches = node_data
  }
  else {
    stop('Value of mutation_position parameter must either "start" or "end"')
  }
  
  apply_res <- apply( no_events_branches, 1,
                      function(elm){
                        tbeta1 <- as.numeric(elm["t_branch_end"])
                        tbeta0 <- as.numeric(elm["t_branch_start"])
                        tbeta1-tbeta0
                      })
  
  x <- sum (node_data$event)/sum( apply_res )
  x
}


# maximum likelihood estimation for lambda (weibull distribution)

lambda_derivative_weib <- function (parms) {
  node_data <- parms$node_data
  p <- parms$p
  mutation_position <- parms$mutation_position
  if (mutation_position == "start"){
    no_events_branches = node_data[node_data["event_indicator"]==0,]
  }
  else if (mutation_position == "end"){
    no_events_branches = node_data
  }
  else {
    stop('Value of mutation_position parameter must either "start" or "end"')
  }
  
  apply_res <- apply( no_events_branches, 1,
                      function(elm){
                        tbeta1 <- as.numeric(elm["t_branch_end"])
                        tbeta0 <- as.numeric(elm["t_branch_start"])
                        if (tbeta1 == 0){0}
                        else {
                          (tbeta1^p)*(1-(tbeta0/tbeta1)^p)
                        }
                      })
  
  x <- (sum (node_data$event)/sum( apply_res ))^(1/p)
  x
  
}


# weibull equation (2) multiplied by p
# root is the maximum likelihood estimation for p

p_derivative <- function (x, parms, draw=FALSE) {
  
  node_data <- parms$node_data
  mutation_position <- parms$mutation_position
  if (mutation_position == "start"){
    no_events_branches = node_data[node_data["event_indicator"]==0,]
  }
  else if (mutation_position == "end"){
    no_events_branches = node_data
  }
  else {
    stop('Value of mutation_position parameter must be either "start" or "end"')
  }
  
  D <- sum (node_data["event_indicator"])
  
  apply_res1 <- apply( no_events_branches, 1,
                       function(elm){ 
                         if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                           tbeta1 <- as.numeric(elm["t_branch_end"])
                           tbeta0 <- as.numeric(elm["t_branch_start"])
                           if (tbeta0 == 0){logtbeta0 <- 0} # mutation can't die before it's birth
                           else {logtbeta0 <- log(tbeta0)}
                           if (tbeta1 == 0 || tbeta1 == 1){ 0 } # from equation 
                           
                           else {
                             logtbeta1 <- log(tbeta1)
                             logtbeta1*(tbeta1^x)*(1-(logtbeta0/logtbeta1)*((tbeta0/tbeta1)^x))
                           }
                         }  
                         else { 0 }
                       })
  
  apply_res2 <- apply( no_events_branches, 1,
                       function(elm){
                         if(!is.na(elm["event_indicator"])){
                           tbeta1 <- as.numeric(elm["t_branch_end"])
                           tbeta0 <- as.numeric(elm["t_branch_start"])
                           if (tbeta1 == 0) { 0 }
                           else {  
                             (tbeta1^x)*(1-(tbeta0/tbeta1)^x)
                           }
                         }
                         else {0}
                       })
  
  apply_res3 <- apply( node_data[node_data["event_indicator"]==1,], 1,
                       function(elm){
                         if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), apply still tries to apply the function to.. header? the list is not empty: it contains the header
                           if (mutation_position == "start"){talpha <- as.numeric(elm["t_branch_start"])}
                           else if (mutation_position == "end"){talpha <- as.numeric(elm["t_branch_end"])}
                           else {stop('Value of mutation_position parameter must be either "start" or "end"')}
                           
                           if (talpha == 0){logtalpha = 0}
                           else {logtalpha = log(talpha)}
                           logtalpha
                         }
                         else {0}
                       }
  )
  
  
  f2 <- D*x*sum( apply_res1 )/sum( apply_res2 ) - x*sum( apply_res3 ) - D
  
  if(draw){
    f2
  }
  else {
    c(F1 = f2)
  }
  
}

# jacobian (just a derivative) of p_derivative_mult function  

p_derivative_jacfunc <- function(x, parms){
  node_data <- parms$node_data
  mutation_position <- parms$mutation_position
  if (mutation_position == "start"){
    no_events_branches = node_data[node_data["event_indicator"]==0,]
  }
  else if (mutation_position == "end"){
    no_events_branches = node_data
  }
  else {
    stop('Value of mutation_position parameter must be either "start" or "end"')
  }
  
  D <- sum (node_data["event_indicator"])
  
  u <- <- apply( no_events_branches, 1,
                 function(elm){ 
                   if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                     tbeta1 <- as.numeric(elm["t_branch_end"])
                     tbeta0 <- as.numeric(elm["t_branch_start"])
                     if (tbeta0 == 0){logtbeta0 <- 0} # mutation can't die before it's birth
                     else {logtbeta0 <- log(tbeta0)}
                     if (tbeta1 == 0 || tbeta1 == 1){ 0 } # from equation 
                     
                     else {
                       logtbeta1 <- log(tbeta1)
                       logtbeta1*(tbeta1^x)*(1-(logtbeta0/logtbeta1)*((tbeta0/tbeta1)^x))
                     }
                   }  
                   else { 0 }
                 })
  
  udot <- apply( no_events_branches, 1,
                 function(elm){ 
                   if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), <apply> still tries to apply the function to.. header? the list is not empty: it contains the header
                     tbeta1 <- as.numeric(elm["t_branch_end"])
                     tbeta0 <- as.numeric(elm["t_branch_start"])
                     if (tbeta0 == 0){logtbeta0 <- 0} # mutation can't die before it's birth
                     else {logtbeta0 <- log(tbeta0)}
                     if (tbeta1 == 0 || tbeta1 == 1){ 0 } # from equation 
                     
                     else {
                       logtbeta1 <- log(tbeta1)
                       (logtbeta1^2)*(tbeta1^x)*(1-(logtbeta0/logtbeta1)^2*((tbeta0/tbeta1)^x))
                     }
                   }  
                   else { 0 }
                 })
  
  v <- apply( no_events_branches, 1,
              function(elm){
                if(!is.na(elm["event_indicator"])){
                  tbeta1 <- as.numeric(elm["t_branch_end"])
                  tbeta0 <- as.numeric(elm["t_branch_start"])
                  if (tbeta1 == 0) { 0 }
                  else {  
                    (tbeta1^x)*(1-(tbeta0/tbeta1)^x)
                  }
                }
                else {0}
              })
  vdot <- apply( no_events_branches, 1,
                 function(elm){
                   if(!is.na(elm["event_indicator"])){
                     tbeta1 <- as.numeric(elm["t_branch_end"])
                     tbeta0 <- as.numeric(elm["t_branch_start"])
                     if (tbeta0 == 0){logtbeta0 <- 0} # mutation can't die before it's birth
                     else {logtbeta0 <- log(tbeta0)}
                     if (tbeta1 == 0 || tbeta1 == 1){ 0 } # from equation
                     else { 
                       logtbeta1 <- log(tbeta1)
                       logtbeta1*(tbeta1^x)*(1-(logtbeta0/logtbeta1)*(tbeta0/tbeta1)^x)
                     }
                   }
                   else {0}
                 })
  
  alphas <- apply( node_data[node_data["event_indicator"]==1,], 1,
                   function(elm){
                     if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), apply still tries to apply the function to.. header? the list is not empty: it contains the header
                       if (mutation_position == "start"){talpha <- as.numeric(elm["t_branch_start"])}
                       else if (mutation_position == "end"){talpha <- as.numeric(elm["t_branch_end"])}
                       else {stop('Value of mutation_position parameter must be either "start" or "end"')}
                       
                       if (talpha == 0){logtalpha = 0}
                       else {logtalpha = log(talpha)}
                       logtalpha
                     }
                     else {0}
                   }
  )
  
  
  fdot <- D*x*(sum(udot)*sum(v) - sum(u)*sum(vdot))/(sum(v))^2 + D*sum(u)/sum(v)- sum(alphas) 
  
}



draw_p_derivative <- function(splitted, anc_node){
  node_data <- splitted[[anc_node]]
  parms <- list(node_data = node_data, mutation_position = "end")
  y <- sapply(seq(from = -5, to = 5, by = 0.05), function (elm){p_derivative(elm,parms, draw=TRUE)})
  plot( x = seq(from = -5, to = 5, by = 0.05), y, type = 'l', xlab = "p", ylab = "f2", axes=F, xaxt="n", yaxt="n", main = anc_node, ylim = c(-10, 5))
  axis(1, pos=0)
  axis(2, pos=0)
  abline(v=0, h=0)
}
