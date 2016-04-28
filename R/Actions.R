prot <- "h1"
prot_data <-  read.csv(paste(c(getwd(), "/input/" ,prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
params <-parameters(splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
##
draw_lnlikelihood (data=splitted, nodename = "169.INTNODE2065", to = 20, by = 0.01, mutation_position = "middle", fishy = TRUE)
draw_hazard(data=splitted, nodename = "78.INTNODE4232", to = 20, by = 0.05, mutation_position = "middle", fishy = TRUE)
##
lrt_all(mutation_position = "middle", fishy = TRUE, tag = "middle_search", pack = "rootsolve", verbose = TRUE)
##
benchmark(parameters(splitted,  jack = FALSE, pack = "nleqslv", filter= FALSE), parameters(splitted,  jack = FALSE, pack = "rootsolve", filter= FALSE),  replications = 1)

##
sink (paste(c(getwd(), "/output/h1_emtest_1clus"),collapse= ""))
em_results <- em_procedure(data=fdata, params=fparams, model = "weibull",iter = 1000, cluster.number= 1, init_method = "cluster", mutation_position = "middle",  filtering = "single", trace = TRUE)
sink()  

sink (paste(c(getwd(), "/output/h1_emtest3"),collapse= ""))
em_results3 <- em_procedure(data=fdata,  params=fparams, iter = 1000, cluster.number= 5, init_method = "by", init_params = c(0.002,8,0.0001,0.9,0.005,0.5,0.002,1.5, 0.0001, 9), init_weights = c(0.15,0.25,0.25, 0.25, 0.1),mutation_position = "middle",  filtering = "single", trace = TRUE)
sink() 




em_results <- em_procedure(data=fdata, params=fparams, iter = 1000, cluster.number= 4, init_method = "by", init_params = c(0.002,8,0.0001,0.7,0.005,0.7,0.002,1.2), init_weights = c(0.25,0.25,0.25, 0.25), mutation_position = "middle",  filtering = "single", trace = TRUE)
