\name{p_numerator_vector}
\alias{p_numerator_vector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p_numerator_vector(no_events_branches, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{no_events_branches}{
%%     ~~Describe \code{no_events_branches} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (no_events_branches, x) 
{
    apply_res <- apply(no_events_branches, 1, function(elm) {
        if (!is.na(elm["event_indicator"])) {
            tbeta1 <- as.numeric(elm["t_branch_end"])
            tbeta0 <- as.numeric(elm["t_branch_start"])
            if (tbeta0 == 0) {
                logtbeta0 <- 0
            }
            else {
                logtbeta0 <- log(tbeta0)
            }
            if (tbeta1 == 0 || tbeta1 == 1) {
                0
            }
            else {
                logtbeta1 <- log(tbeta1)
                logtbeta1 * (tbeta1^x) * (1 - (logtbeta0/logtbeta1) * 
                  ((tbeta0/tbeta1)^x))
            }
        }
        else {
            0
        }
    })
    as.vector(apply_res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
