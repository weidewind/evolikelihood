install.packages("plotly")
## equations for computing maximum likelihood estimators


# maximum likelihood estimation for lambda (exponential distribution)

lambda_derivative_exp <- function(parms){
  data <- parms$data
  rkvector <- parms$rkvector
  if(is.null(rkvector)){
    rkvector <- c(1)
    names(rkvector) = names(data)
  }
  mutation_position <- parms$mutation_position
  
  
  D_vector <- compute_d_vector(data, rkvector)
  D <- sum(D_vector)
  
  branch_vector <- sapply(names(data), function (name){
    with(supplier(data=data,  rkvector=rkvector, elm = name, mutation_position = mutation_position), {
      apply_res <- apply( no_events_branches, 1,
                          function(br){
                            tbeta1 <- as.numeric(br["t_branch_end"])
                            tbeta0 <- as.numeric(br["t_branch_start"])
                            tbeta1-tbeta0
                          })
      apply_res*r
    })
  })

  
  x <- D/sum( branch_vector )
  x
}


# maximum likelihood estimation for lambda (weibull distribution)

lambda_derivative_weib <- function (parms) {
  
  data <- parms$data
  rkvector <- parms$rkvector
  if(is.null(rkvector)){
    rkvector <- c(1)
    names(rkvector) = names(data)
  }
  p <- parms$p
  mutation_position <- parms$mutation_position
  

  
  D_vector <- compute_d_vector(data, rkvector)
  D <- sum(D_vector)
  
  denominator_vector <- sapply(names(data), function (name){
    with(supplier(data=data,  rkvector=rkvector, elm = name, mutation_position = mutation_position), {
      apply_res <- p_denominator_vector(no_events_branches, p)
      apply_res*r
    })
  })
  
  x <- (D/sum(denominator_vector))^(1/p)
  x
  
}

supplier <- function (data,  rkvector, elm, mutation_position){
  node_data <- data[[elm]]
  r <- rkvector[elm]
  if (mutation_position != "start" && mutation_position != "end" && mutation_position != "middle") {
    stop('Value of mutation_position parameter must be either "start" or "end" or "middle"')
  }
  if (mutation_position == "start"){
    no_events_branches = node_data[node_data["event_indicator"]==0,]
  }
  else if (mutation_position == "end"){
    no_events_branches = node_data
  }
  else if (mutation_position == "middle"){
    no_events_branches = node_data[node_data["event_indicator"]==0,]
    first_halves = node_data[node_data["event_indicator"]==1,]
    first_halves[,"t_branch_end"] <- ( first_halves[,"t_branch_end"]+ first_halves[,"t_branch_start"])/2
    no_events_branches = rbind(no_events_branches, first_halves)
  }
  if (mutation_position == "middle"){
    second_halves = node_data[node_data["event_indicator"]==1,]
    second_halves[,"t_branch_start"] <- ( second_halves[,"t_branch_end"]+ second_halves[,"t_branch_start"])/2
    events_branches = second_halves
  }
  else {
    events_branches = node_data[node_data["event_indicator"]==1,]
  }
  
  list(r = r, no_events_branches = no_events_branches, events_branches = events_branches)
}



# weibull equation (2) multiplied by p
# Root is the maximum likelihood estimation for p
# Can be used as model function for multiroot (draw=FALSE)
# or for plotting derivative of loglikelihood function d(logL)/d(p) (weibull equation (2)) (draw=TRUE)

p_derivative <- function (x, parms, draw=FALSE) {
  #data is always a list of dataframes (node_datas)

  mutation_position <- parms$mutation_position
  rkvector <- parms$rkvector
  data <- parms$data
  if(is.null(rkvector)){
    rkvector <- c(1)
    names(rkvector) = names(data)
  }
  


 D_vector <- compute_d_vector(data, rkvector)
 D <- sum(D_vector)
 
 
 numerator_vector <- sapply(names(data), function (name){
   with(supplier(data=data,  rkvector=rkvector, elm = name, mutation_position = mutation_position), {
      apply_res <- p_numerator_vector(no_events_branches, x)
      apply_res*r
  })
 })

 denominator_vector <- sapply(names(data), function (name){
   with(supplier(data=data,  rkvector=rkvector, elm = name, mutation_position = mutation_position), {
     apply_res <- p_denominator_vector(no_events_branches, x)
     apply_res*r
   })
 })

 talpha_vector <- sapply(names(data), function (name){
   with(supplier(data=data,  rkvector=rkvector, elm = name, mutation_position = mutation_position), {
     apply_res <- p_talpha_vector(events_branches, mutation_position)
     apply_res*r
   })
 })

   
 f2 <- D*x*sum(numerator_vector)/sum(denominator_vector) - x*sum(talpha_vector) - D
 #f2 <- Dk*sum(numerator_vector)/sum(denominator_vector) - sum(talpha_vector) - Dk/x
  print (f2)
  if(draw){
    f2
  }
  else {
    c(F1 = f2)
  }

}


# series of subfunctions for computing p derivative (first and second)
compute_d_vector <- function (data, rkvector){
  D_vector <- sapply(names(data), function (elm){
    node_data <- data[[elm]]
    r <- rkvector[elm]
    sum (node_data["event_indicator"]) * r
  })
  D_vector
}


p_numerator_vector <- function(no_events_branches, x){
  apply_res <- apply( no_events_branches, 1,
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
  apply_res
}

p_denominator_vector <-function(no_events_branches, x){
  apply_res <- apply( no_events_branches, 1,
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
  apply_res 
}

p_talpha_vector <- function(events_branches, mutation_position){
  apply_res <- apply( events_branches, 1,
                       function(elm){
                         if(!is.na(elm["event_indicator"])){ #if (nrow(dataframe) == 0), apply still tries to apply the function to.. header? the list is not empty: it contains the header
                           if (mutation_position == "start"){talpha <- as.numeric(elm["t_branch_start"])}
                           else if (mutation_position == "end"){talpha <- as.numeric(elm["t_branch_end"])}
                           else if (mutation_position == "middle"){talpha <- as.numeric(elm["t_branch_start"])} #since events_branches already have start = middle of the branch
                           else {stop('Value of mutation_position parameter must be either "start" or "end" or "middle"')}
                           
                           if (talpha == 0){logtalpha = 0}
                           else {logtalpha = log(talpha)}
                           logtalpha
                         }
                         else {0}
                       })
  apply_res
}

##


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
  events_branches = node_data[node_data["event_indicator"]==1,]
  
  D <- sum (node_data["event_indicator"])
  
  u <- p_numerator_vector(no_events_branches, x)
  
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
  
  v <- p_denominator_vector(no_events_branches, x)
  
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
  
  
  alphas <- p_talpha_vector(events_branches, mutation_position)
  
  
  fdot <- D*x*(sum(udot)*sum(v) - sum(u)*sum(vdot))/(sum(v))^2 + D*sum(u)/sum(v)- sum(alphas) 
  matrix(fdot)
}


# plots partial derivative of loglikelihood function d(logL)/d(p) (weibull equation (2))  for given node_data 

draw_p_derivative <- function(data, mutation_position = "end"){
  #node_data <- splitted[[anc_node]]
  anc_node <- names(data)
  parms <- list(data = data, mutation_position = mutation_position)
  y <- sapply(seq(from = -5, to = 300, by = 0.5), function (elm){p_derivative(elm,parms, draw=TRUE)})
  plot( x = seq(from = -5, to = 300, by = 0.5), y, type = 'l', xlab = "p", ylab = "f2", axes=F, xaxt="n", yaxt="n", main = paste(c(anc_node, " mut pos ", mutation_position)), ylim = c(-10, 5))
  axis(1, pos=0)
  axis(2, pos=0)
  abline(v=0, h=0)
}

draw_lambda_derivative <- function(data, mutation_position = "end"){
  #node_data <- splitted[[anc_node]]
  #node_data <- splitted["138.INTNODE2416"]
  anc_node <- names(data)
 
  y <- sapply(seq(from = -5, to = 20, by = 0.05), function (elm){
    parms <- list(data = data, mutation_position = mutation_position, p=elm)
    lambda_derivative_weib(parms)})
  plot( x = seq(from = -5, to = 20, by = 0.05), y, type = 'l', xlab = "p", ylab = "lambda", axes=F, xaxt="n", yaxt="n", main = anc_node, ylim = c(-0.01, 0.01), xlim = c(-1, 20))
  axis(1, pos=0)
  axis(2, pos=0)
  abline(v=0, h=0)
}

name <-"36.INTNODE1224"
node_data <- splitted[name]
fishy <- TRUE
mutation_position = "middle"

p <- seq(from = 0.025, to = 250, by = 1)
#lambda <- sapply(p, function (elm){
#  parms <- list(data = node_data, mutation_position = mutation_position, p=elm)
#  lambda_derivative_weib(parms)})
lnlikelihood <- sapply(p, function (elm){
      parms <- list(data = node_data, mutation_position = mutation_position, p=elm)
      lambda <- lambda_derivative_weib(parms)
      lnl <- lnlikelihood_weibull(node_data, lambda, elm, fishy = fishy)
      lnl["lnL"]
  }
  )
#surf <- cbind(p, lambda, lnlikelihood )
#plot_ly(x = p, y = lambda, z = lnlikelihood, type = "scatter3D")
#scatterplot3d(p, lambda, lnlikelihood)
plot(p, lnlikelihood, main = paste(c(name, " fishy ", fishy, " mut pos ", mutation_position)))
draw_p_derivative(node_data,  mutation_position= mutation_position)


