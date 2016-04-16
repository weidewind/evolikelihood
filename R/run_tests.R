library('RUnit')
 
source('MLE_functions.R')
 
test.suite <- defineTestSuite("basic",
                              dirs = file.path("tests"),
                              testFileRegexp = '^\\w+\\.R')
 
test.result <- runTestSuite(test.suite)
 
printTextProtocol(test.result)