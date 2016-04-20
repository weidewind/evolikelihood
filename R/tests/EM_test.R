test_initialize_by_clustering <- function (){
  cluster.number = 3
  
  init_values <- initialize_by_clustering(data = splitted, parameters = parameters, cluster.number = cluster.number)
  checkEquals(class(init_values) == "list")
  
  checkEquals(class(init_values$iparameters) == "matrix")
  checkEquals(nrow(init_values$iparameters) == cluster.number)
  checkEquals(ncol(init_values$iparameters) == 2)
  
  checkEquals(length(init_values$iweights) == cluster.number)
  checkEquals(length(init_values$irkvector) == cluster.number)
  
  checkTrue(all(!is.na(init_values$iparameters)))
  checkTrue(all(!is.na(init_values$iweights)))
  checkTrue(all(!is.na(init_values$irkvector)))
  
  list(iparameters = iparameters, iweights = iweights, irkvector = irkvector)
}