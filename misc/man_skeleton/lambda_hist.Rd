\name{lambda_hist}
\alias{lambda_hist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lambda_hist(prot, params, threshold, right = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prot}{
%%     ~~Describe \code{prot} here~~
}
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{right}{
%%     ~~Describe \code{right} here~~
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
function (prot, params, threshold, right = TRUE) 
{
    params <- params[!is.na(params["p_root"]), ]
    p_roots <- params[, "p_root"]
    if (right) {
        lambdas <- params[params["p_root"] > threshold, "lambda_weib_root"]
        sign = ">"
    }
    else {
        lambdas <- params[params["p_root"] <= threshold, "lambda_weib_root"]
        sign = "<="
    }
    h <- hist(lambdas, breaks = 30, plot = FALSE)
    plot(h, main = paste("Histogram of ", prot, " lambda for p", 
        sign, " ", threshold), xlab = "lambda")
    h <- hist(log(lambdas), breaks = 30, plot = FALSE)
    plot(h, main = paste("Histogram of ", prot, " log(lambda) for p", 
        sign, " ", threshold), xlab = "lambda")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
