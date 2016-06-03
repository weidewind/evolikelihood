prot <- "h3"
categories <- 3
model <- "weibull" # weibull or exponential 
em_output <- parse_em_output(prot, model, categories)
## prints parameters and bic for every em launch
em_output
length(em_output)
## 3d picture of model parameters, produced by all launches of em procedure 
## prints bics for each launch
## tells which one is the best
show_em_all_results(em_output, model)
