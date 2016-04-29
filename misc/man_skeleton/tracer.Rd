\name{tracer}
\alias{tracer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tracer(parameters, weights, cluster.number, myplot, init = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{cluster.number}{
%%     ~~Describe \code{cluster.number} here~~
}
  \item{myplot}{
%%     ~~Describe \code{myplot} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
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
function (parameters, weights, cluster.number, myplot, init = FALSE) 
{
    colors <- c("red", "blue", "green", "black", "orange", "gray", 
        "violet")
    if (ncol(parameters) == 2) {
        if (init) {
            myplot <- scatterplot3d(parameters[1, "lambda"], 
                parameters[1, "p"], weights[1], color = colors[1], 
                type = "h", xlim = c(0, 0.1), ylim = c(0, 10), 
                zlim = c(0, 1), pch = 19)
            if (cluster.number > 1) {
                for (i in seq(2, cluster.number, 1)) {
                  myplot$points3d(parameters[i, "lambda"], parameters[i, 
                    "p"], weights[i], col = colors[i], pch = 19, 
                    type = "h")
                }
            }
        }
        else {
            for (i in seq(1, cluster.number, 1)) {
                print("weights")
                print(weights)
                myplot <- myplot
                myplot$points3d(parameters[i, "lambda"], parameters[i, 
                  "p"], weights[i], col = colors[i], type = "h")
            }
        }
    }
    else if (ncol(parameters) == 1) {
        if (init) {
            myplot <- plot(parameters[1, "lambda"], weights[1], 
                col = colors[1], xlim = c(0, 0.1), ylim = c(0, 
                  1), xlab = "lambda", ylab = "weight", pch = 19)
            if (cluster.number > 1) {
                for (i in seq(2, cluster.number, 1)) {
                  points(parameters[i, "lambda"], weights[i], 
                    col = colors[i], pch = 19)
                }
            }
        }
        else {
            for (i in seq(1, cluster.number, 1)) {
                points(parameters[i, "lambda"], weights[i], col = colors[i])
            }
        }
    }
    myplot
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
