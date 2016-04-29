\name{p_derivative_jacfunc}
\alias{p_derivative_jacfunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p_derivative_jacfunc(x, parms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
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
function (x, parms) 
{
    node_data <- parms$node_data
    mutation_position <- parms$mutation_position
    if (mutation_position == "start") {
        no_events_branches = node_data[node_data["event_indicator"] == 
            0, ]
    }
    else if (mutation_position == "end") {
        no_events_branches = node_data
    }
    else {
        stop("Value of mutation_position parameter must be either \"start\" or \"end\"")
    }
    events_branches = node_data[node_data["event_indicator"] == 
        1, ]
    D <- sum(node_data["event_indicator"])
    u <- p_numerator_vector(no_events_branches, x)
    udot <- apply(no_events_branches, 1, function(elm) {
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
                (logtbeta1^2) * (tbeta1^x) * (1 - (logtbeta0/logtbeta1)^2 * 
                  ((tbeta0/tbeta1)^x))
            }
        }
        else {
            0
        }
    })
    v <- p_denominator_vector(no_events_branches, x)
    vdot <- apply(no_events_branches, 1, function(elm) {
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
                  (tbeta0/tbeta1)^x)
            }
        }
        else {
            0
        }
    })
    alphas <- p_talpha_vector(events_branches, mutation_position)
    fdot <- D * x * (sum(udot) * sum(v) - sum(u) * sum(vdot))/(sum(v))^2 + 
        D * sum(u)/sum(v) - sum(alphas)
    matrix(fdot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
