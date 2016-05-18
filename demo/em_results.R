prot <- "h3"
categories <- 5
model <- "exponential" # weibull or exponential 
em_output <- parse_em_output(prot, model, categories)
length(em_output)
show_em_all_results(em_output, model)