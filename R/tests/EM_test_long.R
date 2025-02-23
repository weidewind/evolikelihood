if (!(prot == "h1" && !is.null(splitted) && length(splitted) == 539)){
  prot <- "h1"
  inpath <- paste(c(getwd(), "/data/"),collapse= "")
  prot_data <-  read.csv(paste(c(inpath, prot,"_for_LRT.csv"), collapse=""),stringsAsFactors=FALSE)  
  splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), drop=TRUE)
}
if (is.null(params) || nrow(params) != 142){
  params <- parameters(splitted, mutation_position = "middle", filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
}


## too time-consuming
test_initialize_by_clustering <- function (){
  #arrange
  cluster.number = 3
  #act
  init_values <- initialize_by_clustering(data = splitted, params = params, cluster.number = cluster.number)
  #assert
  checkEquals(class(init_values) == "list")

  checkEquals(class(init_values$iparameters) == "matrix")
  checkEquals(nrow(init_values$iparameters) == cluster.number)
  checkEquals(ncol(init_values$iparameters) == 2)

  checkEquals(length(init_values$iweights) == cluster.number)
  checkEquals(sum(init_values$iweights) == 1)

  checkTrue(all(!is.na(init_values$iparameters)))
  checkTrue(all(!is.na(init_values$iweights)))

}