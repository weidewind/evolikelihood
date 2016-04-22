prot <- "h1"
inpath <- paste(c(getwd(), "/input/"),collapse= "")
prot_data <-  read.csv(paste(c(inpath, prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)


test_compute_d_vector <- function(){

  rkvector <- rep(1, length(splitted))
  names(rkvector) <- names(splitted)
  D_vector<- compute_d_vector(splitted, rkvector)
  checkEquals(sum(D_vector), 524)

}


test_single_mode_supplier <- function() {
  #arrange
  nodename <-  "211.INTNODE1224" 
  #act
  sup <- supplier(data=splitted,  rkvector=c("211.INTNODE1224" = 1), elm = nodename, mutation_position = "middle")
  #assert
  checkTrue(sup$r == 1)
  checkEquals(nrow(sup$no_events_branches), nrow(splitted[[nodename]]))
  checkEquals(nrow(sup$events_branches), 1)
  checkEquals(sup$events_branches["t_branch_start"][[1]], 15)
  checkEquals(nrow(sup$no_events_branches) + nrow(sup$events_branches), 5)
  #total subtree length must stay the same
  checkEquals(sum (sup$events_branches["t_branch_end"] - sup$events_branches["t_branch_start"]) + sum (sup$no_events_branches["t_branch_end"] - sup$no_events_branches["t_branch_start"]), 
            sum(splitted[[nodename]]["t_branch_end"] - splitted[[nodename]]["t_branch_start"]))
}

test_group_mode_supplier <- function() {
  #arrange
  rkvector <- rep(1, length(splitted))
  names(rkvector) <- names(splitted)
  #act
  sup <- supplier(data=splitted,  rkvector=rkvector, elm = "211.INTNODE1224", mutation_position = "middle")
  #assert
  checkTrue(sup$r == 1)
  checkEquals(nrow(sup$no_events_branches), nrow(splitted[[nodename]]))
  checkEquals(nrow(sup$events_branches), 1)
  checkEquals(sup$events_branches["t_branch_start"][[1]], 15)
  checkEquals(nrow(sup$no_events_branches) + nrow(sup$events_branches), 5)
  #total subtree length must stay the same
  checkEquals(sum (sup$events_branches["t_branch_end"] - sup$events_branches["t_branch_start"]) + sum (sup$no_events_branches["t_branch_end"] - sup$no_events_branches["t_branch_start"]), 
              sum(splitted[[nodename]]["t_branch_end"] - splitted[[nodename]]["t_branch_start"]))

}

test_group_mode_lambda_derivative_exp <- function() {
  #arrange
  threenodes <- splitted[c("73.INTNODE1254", "203.INTNODE1737", "169.INTNODE1760")]
  rkvector <- c(1,1,0)
  names(rkvector) <- names(threenodes)
  pars <- list(data = threenodes, rkvector = rkvector, mutation_position = "middle")
  sup73 <- supplier(data=threenodes,  rkvector=rkvector, elm = "73.INTNODE1254", mutation_position = "middle")
  sup203 <- supplier(data=threenodes,  rkvector=rkvector, elm = "203.INTNODE1737", mutation_position = "middle")
  #act
  lambda <- lambda_derivative_exp(pars)
  #assert
  checkEquals(lambda, 2/23)
  
  x <- 1
  
  numsum <- sum(p_numerator_vector(sup73$no_events_branches,x)) + sum(p_numerator_vector(sup203$no_events_branches,x))
  checkEquals(numsum, 8*log(8) + 7.5*log(7.5) - 6*log(6) + 7*log(7) + 2*log(2) + 6.5*log(6.5) - 4*log(4))
  denomsum <-sum(p_denominator_vector(sup73$no_events_branches,x)) + sum(p_denominator_vector(sup203$no_events_branches,x))
  checkEquals(denomsum, 1.5 + 1 +2 +6 + 1 +2.5 + 4 + 3 +2)
  alphasum <- sum(p_talpha_vector(sup73$events_branches, "middle")) + sum(p_talpha_vector(sup203$events_branches, "middle"))
  checkEquals(alphasum, log(7.5) + log(6.5))
  d <- sum(compute_d_vector(threenodes, rkvector))
  checkEquals(d,2)
  p <- p_derivative(x, pars, draw=TRUE)
  checkEquals(p, 2*(8*log(8) + 7.5*log(7.5) - 6*log(6) + 7*log(7) + 2*log(2) + 6.5*log(6.5) - 4*log(4))/23 - log(7.5) - log(6.5)-2)
}