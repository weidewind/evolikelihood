prot <- "h3"
categories <- 3
model <- "weibull" # weibull or exponential 
em_output <- parse_em_output(prot, model, categories)
best <- best_output(em_output)
path <- file.path(getwd(), "output", "wood_likelihood", model, prot, names(best), fsep = .Platform$file.sep)
em_file <-  readLines(path)
zeroline <- tail(grep('^\\s+\\[,1\\]', em_file,  value = FALSE, perl = TRUE), 1)
lowpart <- em_file[seq(zeroline+1,length(em_file),1)]
lastline <- head(grep('^\\[1\\]', lowpart,  value = FALSE, perl = TRUE), 1)
rpart <- lowpart[seq(1, lastline-1, 1)]

h3_antigenic_koel = c(161, 171, 172, 174, 175, 205, 209)
group <- h3_antigenic_koel

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

chisq.test(tbl)

best_output_file <- function (em_output){
  bics <- sapply(em_output, function (elm){
    elm$bic
  })
  sorted_bics <-sort(bics)
  good_names <- names(head(sorted_bics, 1))
  em_output<- em_output[good_names]
}







tbl <- matrix(c(10.3, 3,12.7,20), ncol = 2, byrow = TRUE, nrow = 2)
chisq.test(tbl)