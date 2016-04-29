\name{lambda_derivative_weib}
\alias{lambda_derivative_weib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lambda_derivative_weib(parms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parms}{
%%     ~~Describe \code{parms} here~~
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
function (parms) 
{
    data <- parms$data
    rkvector <- parms$rkvector
    if (is.null(rkvector)) {
        rkvector <- c(1)
        names(rkvector) = names(data)
    }
    p <- parms$p
    mutation_position <- parms$mutation_position
    D_vector <- compute_d_vector(data, rkvector)
    D <- sum(D_vector)
    denominator_vector <- sapply(names(data), function(name) {
        with(supplier(data = data, rkvector = rkvector, elm = name, 
            mutation_position = mutation_position), {
            apply_res <- p_denominator_vector(no_events_branches, 
                p)
            sum(apply_res) * r
        })
    })
    x <- (D/sum(denominator_vector))^(1/p)
    x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
