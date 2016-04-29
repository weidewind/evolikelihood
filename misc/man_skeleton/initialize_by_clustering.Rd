\name{initialize_by_clustering}
\alias{initialize_by_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
initialize_by_clustering(data, params, mutation_position = "middle", cluster.number = 4, model = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
}
  \item{cluster.number}{
%%     ~~Describe \code{cluster.number} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (data, params, mutation_position = "middle", cluster.number = 4, 
    model = NULL) 
{
    if (model == "weibull") {
        parclust <- params[, c("lambda_weib_root", "p_root")]
        scaling <- max(parclust[, 2])/max(parclust[, 1])
        parclust_scaled <- data.frame(lambda_weib_root = scaling * 
            parclust[, 1], p_root = parclust[, 2])
    }
    else {
        parclust_scaled <- data.frame(lambda_exp_root = params[, 
            c("lambda_exp_root")])
    }
    clusters <- kmeans(parclust_scaled, cluster.number)
    categories <- seq(from = 1, to = cluster.number, by = 1)
    rkvectors <- sapply(categories, function(k) {
        rk <- sapply(clusters$cluster, function(category) {
            if (category == k) {
                1
            }
            else {
                0
            }
        })
    })
    rownames(rkvectors) <- params$node
    init_params <- compute_params(data = data, rkvectors = rkvectors, 
        model = model, mutation_position = mutation_position)
    init_weights <- compute_weights(rkvectors)
    list(iparameters = init_params, iweights = init_weights)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
