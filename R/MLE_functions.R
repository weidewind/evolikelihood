list.of.packages <- c("scatterplot3d")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

#Don’t use library() or require(). These modify the search path, affecting what functions are available from the global environment. 
#It’s better to use the DESCRIPTION to specify your package’s requirements
#library(scatterplot3d)
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
 # print (rkvector)
 # print (head(data))
  D_vector <- compute_d_vector(data, rkvector)
  D <- sum(D_vector)
  
 # print ("D")
 # print (D)
  
  branch_vector <- sapply(names(data),  function (name){ #changed apply, 1 to sapply
    with(supplier(data=data,  rkvector=rkvector, elm = name, mutation_position = mutation_position), {
     # print ("r")
     # print (r)
     # print (str(no_events_branches))
      apply_res <- apply( no_events_branches, 1,
                          function(br){
                            #print ("br")
                           # print(class(br))
                           # print (br["t_branch_start"][[1]])
                            tbeta1 <- as.numeric(br["t_branch_end"])
                            tbeta0 <- as.numeric(br["t_branch_start"])
                            #print ( tbeta1-tbeta0)
                            tbeta1-tbeta0
                          })
      sum(apply_res)*r
    })
  })
  
 # print ("sum( branch_vector )")
 # print (sum( branch_vector ))
  
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
      sum(apply_res)*r #changed to sum
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
    no_events_branches <- node_data[node_data["event_indicator"]==0,]
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
      sum(apply_res)*r #changed to sum 
  })
 })

 denominator_vector <- sapply(names(data), function (name){
   with(supplier(data=data,  rkvector=rkvector, elm = name, mutation_position = mutation_position), {
     apply_res <- p_denominator_vector(no_events_branches, x)
     sum(apply_res)*r #changed to sum 
   })
 })

 talpha_vector <- sapply(names(data), function (name){
   with(supplier(data=data,  rkvector=rkvector, elm = name, mutation_position = mutation_position), {
     apply_res <- p_talpha_vector(events_branches, mutation_position)
     sum(apply_res)*r #changed to sum 
   })
 })

   
 f2 <- D*x*sum(numerator_vector)/sum(denominator_vector) - x*sum(talpha_vector) - D
# if (is.nan(f2) || is.na(f2)){
#   print ("Achtung!")
#   print (x)
#   print (sum(numerator_vector))
#   print (sum(denominator_vector))
#   print(sum(talpha_vector))
# }
# else {
#   print ("ok")
#   print (x)
#   print (sum(numerator_vector))
#   print (sum(denominator_vector))
#   print(sum(talpha_vector))
#   print (f2)
# }
 #f2 <- Dk*sum(numerator_vector)/sum(denominator_vector) - sum(talpha_vector) - Dk/x
  #print (f2)
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
  apply_res <- apply( no_events_branches, 1, #changed apply, 1 to saplly
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
  as.vector(apply_res) #changed to as.vector()
}

p_denominator_vector <-function(no_events_branches, x){
  apply_res <- apply( no_events_branches, 1, #changed apply, 1 to saplly
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
  as.vector(apply_res ) #changed to as.vector()
}

p_talpha_vector <- function(events_branches, mutation_position){
  apply_res <- apply( events_branches, 1, #changed apply, 1 to saplly
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
  as.vector(apply_res) #changed to as.vector()
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

draw_p_derivative <- function(data, mutation_position = "end", rkvector, to = 20, by = 0.01){
  #node_data <- splitted[[anc_node]]
  anc_node <- names(data)
  parms <- list(data = data, mutation_position = mutation_position, rkvector = rkvector)
  y <- sapply(seq(from = -5, to = to, by = by), function (elm){p_derivative(elm,parms, draw=TRUE)})
  plot( x = seq(from = -5, to = to, by = by), y, type = 'l', xlab = "p", ylab = "f2", axes=F, xaxt="n", yaxt="n", main = paste(c(anc_node, " mut pos ", mutation_position)), ylim = c(-10, 5))
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

#' @export
#draw_lnlikelihood(splitted, "169.INTNODE1233") 165  node  INTNODE1457 ! 165  node  INTNODE1457 h1 238   INTNODE1364 site  169  node  INTNODE2065
draw_lnlikelihood <- function (data, nodename, to = 20, by = 0.01, mutation_position = "middle", fishy = TRUE){
#data <- splitted
  name <-nodename
  node_data <- data[name]

  p <- seq(from = 0, to = to, by = by)
  lnlikelihood <- sapply(p, function (elm){
      parms <- list(data = node_data, mutation_position = mutation_position, p=elm)
      lambda <- lambda_derivative_weib(parms)
      lnl <- lnlikelihood_weibull(node_data, lambda, elm, fishy = fishy)
      lnl["lnL"]
  }
  )

  plot(p, lnlikelihood, ylab = "lnlikelihood", xlab ="p",  type = 'l', main = paste(c(name, " fishy ", fishy, " mut pos ", mutation_position)))
  draw_p_derivative(node_data,  mutation_position= mutation_position, to=to, by=by)
}

#' @export
draw_hazard <-function(data, nodename, to = 20, by = 0.01, mutation_position = "middle", fishy = TRUE){
  data=splitted
  nodename = "78.INTNODE4232"
  name <-nodename
  node_data <- data[name]
  
  p <- seq(from = 0, to = 30, by = 0.5)
  t <- seq(from = 0, to = 400, by = 5)
 # lambda <- sapply(p, function (elm){
  #  parms <- list(data = node_data, mutation_position = mutation_position, p=elm)
 #   lambda <- lambda_derivative_weib(parms)
 # }
 # )
  hazard <- sapply(t, function (time){
    sapply(p, function (elm){
      parms <- list(data = node_data, mutation_position = mutation_position, p=elm)
      lambda <- lambda_derivative_weib(parms)
      elm*(lambda^elm)*(time^(elm-1))
    })
  })
 hazard_time <-sapply(t, function (time){

     parms <- list(data = node_data, mutation_position = mutation_position, p=30)
     lambda <- lambda_derivative_weib(parms)
     30*(lambda^30)*(time^(30-1))

 })
 # scatterplot3d()
  persp(p, t, hazard, zlim=c(0,1),phi = 45, theta = 45,
        xlab = "p", ylab = "t",
        main = "hazard"
  )
 
 plot(t, hazard_time, ylab = "hazard", xlab ="t",  type = 'l', ylim = c(0, 100), main = paste(c(name, " fishy ", fishy, " mut pos ", mutation_position)))
 # draw_p_derivative(node_data,  mutation_position= mutation_position, to=to, by=by)
}

