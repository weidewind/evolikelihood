\name{initialize_by}
\alias{initialize_by}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
initialize_by(init_params, init_weights, model = NULL, cluster.number = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init_params}{
%%     ~~Describe \code{init_params} here~~
}
  \item{init_weights}{
%%     ~~Describe \code{init_weights} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{cluster.number}{
%%     ~~Describe \code{cluster.number} here~~
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
function (init_params, init_weights, model = NULL, cluster.number = 4) 
{
    if (class(init_params) != "numeric") {
        stop("Invalid argument params: expected vector of length 2*cluster.number")
    }
    if (model == "weibull") {
        if (length(init_params) != 2 * cluster.number) {
            stop("Invalid params length: expected vector of length 2*cluster.number")
        }
    }
    else if (model == "exponential") {
        if (length(init_params) != cluster.number) {
            stop("Invalid params length: expected vector of length cluster.number")
        }
    }
    else {
        stop(paste(c("Invalid model name: expected 'weibull' or 'exponential', recieved ", 
            model), collpase = ""))
    }
    if (length(init_weights) != cluster.number) {
        stop("Invalid weights length: expected vector of length cluster.number")
    }
    if (sum(init_weights) != 1) {
        stop("Invalid weights value: sum of weights must be equal to 1")
    }
    categories <- seq(from = 1, to = cluster.number, by = 1)
    if (model == "weibull") {
        init_params <- matrix(init_params, ncol = 2, nrow = cluster.number, 
            byrow = TRUE)
        colnames(init_params) = c("lambda", "p")
    }
    else if (model == "exponential") {
        init_params <- matrix(init_params, ncol = 1, nrow = cluster.number, 
            byrow = TRUE)
        colnames(init_params) = c("lambda")
    }
    list(iparameters = init_params, iweights = init_weights)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
