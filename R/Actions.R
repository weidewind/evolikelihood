prot <- "h1"
prot_data <-  read.csv(paste(c(getwd(), "/input/" ,prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
params <-parameters(splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)

draw_lnlikelihood (data=splitted, nodename = "169.INTNODE2065", to = 20, by = 0.01, mutation_position = "middle", fishy = TRUE)
draw_hazard(data=splitted, nodename = "78.INTNODE4232", to = 20, by = 0.05, mutation_position = "middle", fishy = TRUE)

lrt_all(mutation_position = "middle", fishy = TRUE, tag = "middle_search", pack = "rootsolve", verbose = TRUE)

benchmark(parameters(splitted,  jack = FALSE, pack = "nleqslv", filter= FALSE), parameters(splitted,  jack = FALSE, pack = "rootsolve", filter= FALSE),  replications = 1)
