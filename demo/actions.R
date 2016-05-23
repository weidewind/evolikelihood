list.of.packages <- c("rbenchmark")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')


prot <- "h1"
prot_data <-  read.csv(paste(c(getwd(), "/data/" ,prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
params <-parameters(splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)

proth3 <- "h3"
prot_datah3 <-  read.csv(paste(c(getwd(), "/data/" ,proth3,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splittedh3 <- split(prot_datah3, list(prot_datah3$site, prot_datah3$ancestor_node), drop=TRUE)
paramsh3 <-parameters(splittedh3, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)

protn1 <- "n1"
prot_datan1 <-  read.csv(paste(c(getwd(), "/data/" ,protn1,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splittedn1 <- split(prot_datan1, list(prot_datan1$site, prot_datan1$ancestor_node), drop=TRUE)
paramsn1 <-parameters(splittedn1, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)

protn2 <- "n2"
prot_datan2 <-  read.csv(paste(c(getwd(), "/data/" ,protn2,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splittedn2 <- split(prot_datan2, list(prot_datan2$site, prot_datan2$ancestor_node), drop=TRUE)
paramsn2 <-parameters(splittedn2, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)






#params <-parameters(splitted, mutation_position = "middle",  filter = FALSE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
#dumb_wood_likelihood(data = splitted, prot = prot, mutation_position = "middle", tag = "all_subtrees", fishy=TRUE, params = params, threshold = 0, all = TRUE)
dumb_wood_likelihood(data = splittedh3, prot = proth3, mutation_position = "middle", tag = "all_subtrees", fishy=TRUE, params = paramsh3, threshold = 0, all = TRUE)
dumb_wood_likelihood(data = splittedn2, prot = protn2, mutation_position = "middle", tag = "all_subtrees", fishy=TRUE, params = paramsn2, threshold = 0, all = TRUE)
dumb_wood_likelihood(data = splittedn1, prot = protn1, mutation_position = "middle", tag = "all_subtrees", fishy=TRUE, params = paramsn1, threshold = 0, all = TRUE)




##
draw_lnlikelihood (data=splitted, nodename = "169.INTNODE2065", to = 20, by = 0.01, mutation_position = "middle", fishy = TRUE)
draw_hazard(data=splitted, nodename = "78.INTNODE4232", to = 20, by = 0.05, mutation_position = "middle", fishy = TRUE)
##


##
lrt_all(mutation_position = "middle", fishy = TRUE, tag = "all_subtrees", pack = "rootsolve", verbose = TRUE, threshold = 0, all = TRUE)
##
benchmark(parameters(splitted,  jack = FALSE, pack = "nleqslv", filter= FALSE), parameters(splitted,  jack = FALSE, pack = "rootsolve", filter= FALSE),  replications = 1)

##
sink (paste(c(getwd(), "/output/h1_emtest_2cl"),collapse= ""))
em_results <- em_procedure(data=splitted, params=params, model = "weibull",iter = 1000, cluster.number= 2, init_method = "cluster", mutation_position = "middle",  filtering = "single", trace = TRUE)
sink()  

sink (paste(c(getwd(), "/output/h1_emtest3"),collapse= ""))
em_results3 <- em_procedure(data=fdata,  params=fparams, iter = 1000, cluster.number= 5, init_method = "by", init_params = c(0.002,8,0.0001,0.9,0.005,0.5,0.002,1.5, 0.0001, 9), init_weights = c(0.15,0.25,0.25, 0.25, 0.1),mutation_position = "middle",  filtering = "single", trace = TRUE)
sink() 




em_results <- em_procedure(data=fdata, params=fparams, iter = 1000, cluster.number= 4, init_method = "by", init_params = c(0.002,8,0.0001,0.7,0.005,0.7,0.002,1.2), init_weights = c(0.25,0.25,0.25, 0.25), mutation_position = "middle",  filtering = "single", trace = TRUE)
