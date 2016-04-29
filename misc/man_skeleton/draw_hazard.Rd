\name{draw_hazard}
\alias{draw_hazard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
draw_hazard(data, nodename, to = 20, by = 0.01, mutation_position = "middle", fishy = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{nodename}{
%%     ~~Describe \code{nodename} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
}
  \item{by}{
%%     ~~Describe \code{by} here~~
}
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
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
function (data, nodename, to = 20, by = 0.01, mutation_position = "middle", 
    fishy = TRUE) 
{
    data = splitted
    nodename = "78.INTNODE4232"
    name <- nodename
    node_data <- data[name]
    p <- seq(from = 0, to = 30, by = 0.5)
    t <- seq(from = 0, to = 400, by = 5)
    hazard <- sapply(t, function(time) {
        sapply(p, function(elm) {
            parms <- list(data = node_data, mutation_position = mutation_position, 
                p = elm)
            lambda <- lambda_derivative_weib(parms)
            elm * (lambda^elm) * (time^(elm - 1))
        })
    })
    hazard_time <- sapply(t, function(time) {
        parms <- list(data = node_data, mutation_position = mutation_position, 
            p = 30)
        lambda <- lambda_derivative_weib(parms)
        30 * (lambda^30) * (time^(30 - 1))
    })
    persp(p, t, hazard, zlim = c(0, 1), phi = 45, theta = 45, 
        xlab = "p", ylab = "t", main = "hazard")
    plot(t, hazard_time, ylab = "hazard", xlab = "t", type = "l", 
        ylim = c(0, 100), main = paste(c(name, " fishy ", fishy, 
            " mut pos ", mutation_position)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
