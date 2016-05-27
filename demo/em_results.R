prot <- "n1"
categories <- 1
model <- "weibull" # weibull or exponential 
#em_output <- parse_em_output(prot, model, categories)
em_output <- parse_group_output(prot, model, categories)
em_output
length(em_output)
show_em_all_results(em_output, model)