\name{draw_p_derivative}
\alias{draw_p_derivative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
draw_p_derivative(data, mutation_position = "end", rkvector, to = 20, by = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
}
  \item{rkvector}{
%%     ~~Describe \code{rkvector} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
}
  \item{by}{
%%     ~~Describe \code{by} here~~
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
function (data, mutation_position = "end", rkvector, to = 20, 
    by = 0.01) 
{
    anc_node <- names(data)
    parms <- list(data = data, mutation_position = mutation_position, 
        rkvector = rkvector)
    y <- sapply(seq(from = -5, to = to, by = by), function(elm) {
        p_derivative(elm, parms, draw = TRUE)
    })
    plot(x = seq(from = -5, to = to, by = by), y, type = "l", 
        xlab = "p", ylab = "f2", axes = F, xaxt = "n", yaxt = "n", 
        main = paste(c(anc_node, " mut pos ", mutation_position)), 
        ylim = c(-10, 5))
    axis(1, pos = 0)
    axis(2, pos = 0)
    abline(v = 0, h = 0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
