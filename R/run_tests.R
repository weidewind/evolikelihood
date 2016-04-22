library('RUnit')

dirpath <- paste(c(getwd(), "/R/"),collapse= "")
source(paste(c(dirpath, "MLE_functions.R"),collapse= ""))
source(paste(c(dirpath, "MLE_parameters.R"),collapse= ""))
source(paste(c(dirpath, "LRT.R"),collapse= ""))
source(paste(c(dirpath, "EM.R"),collapse= ""))


test.suite <- defineTestSuite("basic",
                              #dirs = file.path("tests"),
                              paste(c(dirpath, "tests"),collapse= ""),
                              testFileRegexp = 's_.R')
 
test.result <- runTestSuite(test.suite)
 
printTextProtocol(test.result)