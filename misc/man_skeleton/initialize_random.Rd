\name{initialize_random}
\alias{initialize_random}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
initialize_random(params, cluster.number = 4, model = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{cluster.number}{
%%     ~~Describe \code{cluster.number} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (params, cluster.number = 4, model = NULL) 
{
    if (model == "weibull") {
        init_params <- matrix(ncol = 2, nrow = cluster.number, 
            byrow = TRUE)
        colnames(init_params) = c("lambda", "p")
        p_roots <- params[, "p_root"]
        prand <- exp(runif(cluster.number, min = min(log(p_roots)), 
            max = max(log(p_roots))))
        init_params[, "p"] <- prand
        lambda_roots <- params[, "lambda_weib_root"]
        lambdarand <- exp(runif(cluster.number, min = min(log(lambda_roots)), 
            max = max(log(lambda_roots))))
        init_params[, "lambda"] <- lambdarand
    }
    else if (model == "exponential") {
        init_params <- matrix(ncol = 1, nrow = cluster.number, 
            byrow = TRUE)
        colnames(init_params) = c("lambda")
        lambda_roots <- params[, "lambda_exp_root"]
        lambdarand <- exp(runif(cluster.number, min = min(log(lambda_roots)), 
            max = max(log(lambda_roots))))
        init_params[, "lambda"] <- lambdarand
    }
    init_weights = rep(1/cluster.number, cluster.number)
    list(iparameters = init_params, iweights = init_weights)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
