library('RUnit')

dirpath <- paste(c(getwd(), "/R/"),collapse= "")
source(paste(c(dirpath, "MLE_functions.R"),collapse= ""))
source(paste(c(dirpath, "MLE_parameters.R"),collapse= ""))
source(paste(c(dirpath, "LRT.R"),collapse= ""))
source(paste(c(dirpath, "EM.R"),collapse= ""))


test.suite <- defineTestSuite("basic",
                              #dirs = file.path("tests"),
                              paste(c(dirpath, "tests"),collapse= ""),
                              testFileRegexp = 'test.R')
 
test.result <- runTestSuite(test.suite)
 
printTextProtocol(test.result)

test.suite.long <- defineTestSuite("long",
                              #dirs = file.path("tests"),
                              paste(c(dirpath, "tests"),collapse= ""),
                              testFileRegexp = 'test_long.R')

test.result.long <- runTestSuite(test.suite.long)

printTextProtocol(test.result.long)