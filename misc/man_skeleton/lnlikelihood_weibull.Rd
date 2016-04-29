\name{lnlikelihood_weibull}
\alias{lnlikelihood_weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lnlikelihood_weibull(node_data, lambda, p, fishy = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node_data}{
%%     ~~Describe \code{node_data} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{fishy}{
%%     ~~Describe \code{fishy} here~~
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
function (node_data, lambda, p, fishy = FALSE) 
{
    if (class(node_data) == "list") {
        node_data <- node_data[[1]]
    }
    beta_branches = node_data[node_data["event_indicator"] == 
        0, ]
    alpha_branches = node_data[node_data["event_indicator"] == 
        1, ]
    if (fishy) {
        apply_res1 <- apply(alpha_branches, 1, function(elm) {
            if (!is.na(elm["event_indicator"])) {
                talpha1 <- as.numeric(elm["t_branch_end"])
                talpha0 <- as.numeric(elm["t_branch_start"])
                talpha_middle <- (talpha1 + talpha0)/2
                log(p) + p * log(lambda) + (p - 1) * log(talpha_middle) + 
                  log(talpha1 - talpha_middle)
            }
            else {
                0
            }
        })
    }
    else {
        apply_res1 <- apply(alpha_branches, 1, function(elm) {
            if (!is.na(elm["event_indicator"])) {
                talpha1 <- as.numeric(elm["t_branch_end"])
                talpha0 <- as.numeric(elm["t_branch_start"])
                survival <- exp((lambda^p) * (talpha1^p) * (((talpha0/talpha1)^p) - 
                  1))
                log(1 - survival)
            }
            else {
                0
            }
        })
    }
    apply_res2 <- apply(beta_branches, 1, function(elm) {
        if (!is.na(elm["event_indicator"])) {
            tbeta1 <- as.numeric(elm["t_branch_end"])
            tbeta0 <- as.numeric(elm["t_branch_start"])
            if (tbeta1 == 0) {
                0
            }
            else {
                (tbeta1^p) * (((tbeta0/tbeta1)^p) - 1)
            }
        }
        else {
            0
        }
    })
    if (fishy) {
        apply_res3 <- apply(alpha_branches, 1, function(elm) {
            if (!is.na(elm["event_indicator"])) {
                talpha1 <- as.numeric(elm["t_branch_end"])
                talpha0 <- as.numeric(elm["t_branch_start"])
                talpha_middle <- (talpha1 + talpha0)/2
                (talpha_middle^p) * (((talpha0/talpha_middle)^p) - 
                  1)
            }
            else {
                0
            }
        })
    }
    if (fishy) {
        lnL <- sum(apply_res1) + (lambda^p) * (sum(apply_res2) + 
            sum(apply_res3))
    }
    else {
        lnL <- sum(apply_res1) + (lambda^p) * sum(apply_res2)
    }
    c(lnL = lnL, AIC = aic(lnL, 2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
