inpath <- paste(c(getwd(), "/input/"),collapse= "")

prot <- "h3"
prot_data <-  read.csv(paste(c(inpath, prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)



test_supplier <- function() {
  nodename <-  "211.INTNODE1224" 
  sup <- supplier(data=splitted,  rkvector=c("211.INTNODE1224" = 1), elm = nodename, mutation_position = "middle")
  checkTrue(sup$r == 1)
  checkEquals(nrow(sup$no_events_branches), nrow(splitted[[nodename]]))
  checkEquals(nrow(sup$events_branches), 1)
  checkEquals(sup$events_branches["t_branch_start"][[1]], 15)
  checkEquals(nrow(sup$no_events_branches) + nrow(sup$events_branches), 5)
  #total subtree length must stay the same
  checkEquals(sum (sup$events_branches["t_branch_end"] - sup$events_branches["t_branch_start"]) + sum (sup$no_events_branches["t_branch_end"] - sup$no_events_branches["t_branch_start"]), 
            sum(splitted[[nodename]]["t_branch_end"] - splitted[[nodename]]["t_branch_start"]))

#list(r = 1, no_events_branches = no_events_branches, events_branches = events_branches)
   
#isTRUE(all.equal(target, current, ...,
#          check.attributes = TRUE, use.names = TRUE)) # for checking equality of two lists

#   checkEquals(divideBy(4, 2), 2)
#    checkTrue(is.na(divideBy(4, 0)))
#    checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}