\name{lrt}
\alias{lrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lrt(node_data, lambda_exp, lambda_weib, p, fishy = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node_data}{
%%     ~~Describe \code{node_data} here~~
}
  \item{lambda_exp}{
%%     ~~Describe \code{lambda_exp} here~~
}
  \item{lambda_weib}{
%%     ~~Describe \code{lambda_weib} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{fishy}{
%%     ~~Describe \code{fishy} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (node_data, lambda_exp, lambda_weib, p, fishy = FALSE, 
    verbose = FALSE) 
{
    weibull <- lnlikelihood_weibull(node_data, lambda = lambda_weib, 
        p = p, fishy)
    exponential <- lnlikelihood_exp(node_data, lambda = lambda_exp, 
        fishy)
    LR <- 2 * (weibull[1] - exponential[1])
    if (verbose) {
        print("weibull")
        print(weibull)
        print("exponential")
        print(exponential)
        print("LR")
        print(LR)
    }
    LR
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
