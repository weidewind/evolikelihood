\name{compute_model_lnL}
\alias{compute_model_lnL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_model_lnL(data, model = NULL, parameters, weights)
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
function (data, model = NULL, parameters, weights) 
{
    cluster.number = length(weights)
    categories <- seq(from = 1, to = cluster.number, by = 1)
    likelihood_vector <- sapply(names(data), function(nodename) {
        cat_probs <- sapply(categories, function(cat) {
            if (model == "weibull") {
                lnL_dat <- lnlikelihood_weibull(data[[nodename]], 
                  parameters[cat, "lambda"], parameters[cat, 
                    "p"], fishy = TRUE)
            }
            else {
                lnL_dat <- lnlikelihood_exp(data[[nodename]], 
                  parameters[cat, "lambda"], fishy = TRUE)
            }
            lnL <- lnL_dat[1]
            weights[cat] * exp(lnL)
        })
        sum(cat_probs)
    })
    lnL <- sum(log(likelihood_vector))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
