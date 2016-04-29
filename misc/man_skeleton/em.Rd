\name{em}
\alias{em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
em(data, model = NULL, parameters, weights, iter = 100, cluster.number = 4, mutation_position = "middle", trace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{cluster.number}{
%%     ~~Describe \code{cluster.number} here~~
}
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
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
function (data, model = NULL, parameters, weights, iter = 100, 
    cluster.number = 4, mutation_position = "middle", trace = TRUE) 
{
    if (trace) {
        myplot <- tracer(parameters, weights, cluster.number, 
            init = TRUE)
    }
    old_lnL <- NULL
    for (i in seq(1, iter, 1)) {
        print(paste(c("------------Step ", i), collapse = ""))
        rkvectors <- compute_rkvectors(data = data, parameters = parameters, 
            model = model, weights = weights)
        parameters <- compute_params(data = data, rkvectors = rkvectors, 
            model = model, mutation_position = mutation_position)
        weights <- compute_weights(rkvectors)
        print(rkvectors)
        print(weights)
        if (trace) {
            myplot <- tracer(parameters, weights, cluster.number, 
                myplot = myplot, init = FALSE)
        }
        model_lnL <- compute_model_lnL(data = data, model = model, 
            parameters = parameters, weights = weights)
        print("model lnL")
        print(model_lnL)
        if (!is.null(old_lnL) && model_lnL - old_lnL < 1e-04) {
            break
        }
        else {
            old_lnL <- model_lnL
        }
    }
    model_bic <- bic(lnL = model_lnL, model = model, cluster.number = cluster.number, 
        n = length(data))
    print("data length")
    print(length(data))
    print("model bic")
    print(model_bic)
    list(parameters = parameters, rkvectors = rkvectors, weights = weights, 
        lnL = model_lnL, bic = model_bic)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
