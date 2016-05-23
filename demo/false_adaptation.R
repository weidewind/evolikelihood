prot <- "h3"
model <- "exponential"
categories <- 1

prot_data <-  read.csv(file.path(getwd(), "data", paste(c(prot,"_for_LRT.csv"), collapse=""),fsep = .Platform$file.sep),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)

group <- c("236.INTNODE2473", "189.INTNODE3183", "241.INTNODE3625",  "66.INTNODE3936", "218.INTNODE3933")

grepper <- paste(group,collapse="|")
group_nodes <- grep(grepper, names(splitted))
group_splitted <- splitted[group_nodes]
params <-parameters(group_splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
sink (file.path(getwd(), "output","group_likelihood", model, prot, "false_adaptation", fsep = .Platform$file.sep))
em_results <- em_procedure(data=group_splitted, params=params, model = model,  cluster.number= 1, init_method = "cluster", mutation_position = "middle",  filtering = "single", trace = FALSE)
sink() 