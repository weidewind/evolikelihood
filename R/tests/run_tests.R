list.of.packages <- c("RUnit")
new.packages <- setdiff(list.of.packages, installed.packages()[,"Package"])
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library('RUnit')
library(evolike) #instead of sourcing multiple files

dirpath <- paste(c(getwd(), "/R/"),collapse= "")
#source(paste(c(dirpath, "MLE_functions.R"),collapse= ""))
#source(paste(c(dirpath, "MLE_parameters.R"),collapse= ""))
#source(paste(c(dirpath, "LRT.R"),collapse= ""))
#source(paste(c(dirpath, "EM.R"),collapse= ""))


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