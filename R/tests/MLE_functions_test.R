prot <- "h1"
prot_data <-  read.csv(paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/likelihood/nsyn/",prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)



test_divideBy <- function() {
supplier(data=splitted,  rkvector=c(name = 1), elm = name, mutation_position = "middle"), list(r = r, no_events_branches = no_events_branches, events_branches = events_branches)
   
isTRUE(all.equal(target, current, ...,
          check.attributes = TRUE, use.names = TRUE)) -  для сравнения листов. если объекты не одинаковы, all.equal выдаст различия, так что нельзя его без is.TRUE использовать в ифах 

   checkEquals(divideBy(4, 2), 2)
    checkTrue(is.na(divideBy(4, 0)))
    checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}