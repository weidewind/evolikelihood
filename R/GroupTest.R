
# makes em procedure for groups
wood_groups_test <- function(prot, model, tag){
  data("mygroups")
  
  
  #prot <- "h3"
  prot_groups <-mygroups[grep(paste(c("^", prot, "_"), collapse = ""), names(mygroups))]
  #model <- "weibull"
  categories <- 1
  
  prot_data <-  read.csv(file.path(getwd(), "data", paste(c(prot,"_for_LRT.csv"), collapse=""),fsep = .Platform$file.sep),stringsAsFactors=FALSE)  
  splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
  
  count_cores <- detectCores() - 1
  # Initiate cluster
  cl <- makeCluster(count_cores)
  clusterExport(cl, list("prot", "prot_groups", "splitted", "model"), envir = environment())
  clusterCall(cl, function() library(evolike))
  
  prot_groups_names <- names(prot_groups)
  em_results_list <- parLapply(cl, prot_groups_names, function(group_name){
    sink (file.path(getwd(), "output","group_likelihood", model, prot, paste(c("logs_", group_name), collapse= ""), tag, fsep = .Platform$file.sep))
    
    print (group_name)
    group <- prot_groups[group_name][[1]]
    group_pattern <- sapply(group, function (e){
      paste(c("^", e, "\\."), collapse = "")
    })
    grepper <- paste(group_pattern,collapse="|")
    group_nodes <- grep(grepper, names(splitted))
    group_splitted <- splitted[group_nodes]
    
    compl <- !(seq(1,length(splitted),1) %in% group_nodes)
    compl_splitted <- splitted[compl]
    
    params <-parameters(group_splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
    compl_params <-parameters(compl_splitted, mutation_position = "middle",  filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
    
    sink() 
    
    sink (file.path(getwd(), "output","group_likelihood", model, prot, group_name, fsep = .Platform$file.sep))
    em_results <- em_procedure(data=group_splitted, params=params, model = model,  cluster.number= 1, init_method = "cluster", mutation_position = "middle",  filtering = "single", trace = FALSE)
    sink() 
    sink (file.path(getwd(), "output","group_likelihood", model, prot, paste(c(group_name, "_complement"), collapse= ""), fsep = .Platform$file.sep))
    em_results <- em_procedure(data=compl_splitted, params=compl_params, model = model,  cluster.number= 1, init_method = "cluster", mutation_position = "middle",  filtering = "single", trace = FALSE)
    sink() 
    
    
    em_results
  })
  
  
  
  
  stopCluster(cl)
}


## takes a posterior probabilities from wood_likelihood analysis (you must explicitly select the optimal number of categories)
chisq_groups_test < -function(prot, categories, tag){
  data("mygroups")
  model <- "weibull"
  prot_groups <-mygroups[grep(paste(c("^", prot, "_"), collapse = ""), names(mygroups))]
  
  em_output <- parse_em_output(prot, model, categories)
  best <- best_output(em_output)
  path <- file.path(getwd(), "output", "wood_likelihood", model, prot, names(best), fsep = .Platform$file.sep)
  em_file <-  readLines(path)
  zeroline <- tail(grep('^\\s+\\[,1\\]', em_file,  value = FALSE, perl = TRUE), 1)
  lowpart <- em_file[seq(zeroline+1,length(em_file),1)]
  lastline <- head(grep('^\\[1\\]', lowpart,  value = FALSE, perl = TRUE), 1)
  rpart <- lowpart[seq(1, lastline-1, 1)]
  
  
  count_cores <- detectCores() - 1
  # Initiate cluster
  cl <- makeCluster(count_cores)
  clusterExport(cl, list("prot", "rpart", "prot_groups", "splitted", "model", "categories", "best"), envir = environment())
  clusterCall(cl, function() library(evolike))
  
  prot_groups_names <- names(prot_groups)
  em_results_list <- parLapply(cl, prot_groups_names, function(group_name){
    sink (file.path(getwd(), "output", "chisq", model, prot, group_name, tag, fsep = .Platform$file.sep))
    
    group <- prot_groups[group_name][[1]]
    rvectors <- lapply(rpart, function(row){
      splittedrow <- unlist(strsplit(row, '\\s+'))
      my_params <- lapply( tail(splittedrow, length(splittedrow)), function(elm){
        elm
      }) 
      
    })
    
    df <- data.frame(matrix(unlist(rvectors), nrow=length(rvectors), byrow=T),stringsAsFactors=FALSE)
    rownames(df) <- df[,1]
    df <- df[,seq(2, ncol(df), 1)]
    
    group_pattern <- sapply(group, function (e){
      paste(c("^", e, "\\."), collapse = "")
    })
    grepper <- paste(group_pattern,collapse="|")
    group_nodes <- grep(grepper, rownames(df))
    group_df <- df[group_nodes,]
    
    compl_nodes <- !(seq(1,nrow(df),1) %in% group_nodes)
    compl_df <- df[compl_nodes,]
    
    ageing <- sapply(seq(1, categories, 1), function(cat){
      if (best[[1]]$params[cat,"p"] > 1){TRUE} else {FALSE}
    })
    
    
    ageing_g <- sum(sapply(group_df[,(ageing)], function(e){sum(as.numeric(e))}))
    adapt_g <- sum(sapply(group_df[,!(ageing)], function(e){sum(as.numeric(e))}))
    ageing_c <- sum(sapply(compl_df[,(ageing)], function(e){sum(as.numeric(e))}))
    adapt_c <- sum(sapply(compl_df[,!(ageing)], function(e){sum(as.numeric(e))}))
    
    tbl <- matrix(c(ageing_g, adapt_g, ageing_c, adapt_c), ncol = 2, byrow = TRUE, nrow = 2)
    colnames(tbl) <- c("ageing", "adapt")
    rownames(tbl) <- c("group", "compl")
    print(tbl)
    
    print(chisq.test(tbl))
    print(fisher.test(tbl))
    sink() 
    
    
  })
  
  stopCluster(cl)

}


best_output_file <- function (em_output){
  bics <- sapply(em_output, function (elm){
    elm$bic
  })
  sorted_bics <-sort(bics)
  good_names <- names(head(sorted_bics, 1))
  em_output<- em_output[good_names]
}


