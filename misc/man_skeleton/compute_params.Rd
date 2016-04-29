\name{compute_params}
\alias{compute_params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_params(data, model = NULL, rkvectors, mutation_position = "middle")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{rkvectors}{
%%     ~~Describe \code{rkvectors} here~~
}
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, model = NULL, rkvectors, mutation_position = "middle") 
{
    cluster.number = ncol(rkvectors)
    categories <- seq(from = 1, to = cluster.number, by = 1)
    if (model == "weibull") {
        new_params <- matrix(nrow = cluster.number, ncol = 2)
        colnames(new_params) = c("lambda", "p")
    }
    else {
        new_params <- matrix(nrow = cluster.number, ncol = 1)
        colnames(new_params) = c("lambda")
    }
    for (k in categories) {
        k_params <- find_single_root(data = data, mutation_position = mutation_position, 
            rkvector = rkvectors[, k], jack = FALSE, pack = "rootsolve", 
            verbose = TRUE)
        if (model == "weibull") {
            new_params[k, "p"] <- k_params["p_root"]
            new_params[k, "lambda"] <- k_params["lambda_weib_root"]
        }
        else {
            new_params[k, "lambda"] <- k_params["lambda_exp_root"]
        }
    }
    new_params
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
