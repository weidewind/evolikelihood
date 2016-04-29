\name{bic}
\alias{bic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bic(lnL, model = NULL, cluster.number, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lnL}{
%%     ~~Describe \code{lnL} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{cluster.number}{
%%     ~~Describe \code{cluster.number} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (lnL, model = NULL, cluster.number, n) 
{
    if (model == "weibull") {
        p <- 3 * cluster.number - 1
    }
    else {
        p <- 2 * cluster.number - 1
    }
    -2 * lnL + p * log(n)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
