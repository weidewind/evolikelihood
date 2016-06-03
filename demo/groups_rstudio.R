prot <- "h3"
model <- "weibull"
categories <- 1
## prints parameters and bic for every group. Never tested for number of categories > 1
em_output <- parse_group_output(prot, model, categories)
em_output

## group likelihoods for exp and weib models, almost the same procedure as for wood likelihood, 
## except that we compute parameters for 1 category only 
## (no optimization, do not be confused by em_procedure that is used here )
## output  - in output/group_likelihood
prot <- "h3"
model <- "weibull" # weibull or exponential 
wood_groups_test (prot, model, "mytag")

## LRT for those data (won't process files with tag)
parse_group_LRT("h3")

## another way to test groups: chisq on contingency table produced from wood_likelihood/weibull data
## compares sums of a posteriori probabilities to belong to categories with p>1 and p <1  - for group and for its complement
## Attention! you must manually select the optimal number of categories for the protein (taken from wood_likelihood.xls, sheet wood_likelihood)
## output  - in output/chisq

prot <- "n2"
categories <- 3
chisq_groups_test (prot, categories, "mytag")

