library('RUnit')

dirpath <- paste(c(getwd(), "/R/"),collapse= "")
source(paste(c(dirpath, "MLE_functions.R"),collapse= ""))
 
test.suite <- defineTestSuite("basic",
                              #dirs = file.path("tests"),
                              paste(c(dirpath, "tests"),collapse= ""),
                              testFileRegexp = '.R')
 
test.result <- runTestSuite(test.suite)
 
printTextProtocol(test.result)