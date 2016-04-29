\name{em_procedure}
\alias{em_procedure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
em_procedure(data, params, model = NULL, iter = 100, cluster.number = 4, init_method = c("cluster", "random", "by"), init_params = NULL, init_weights = NULL, mutation_position = "middle", filtering = c("single", "unsolved"), trace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{cluster.number}{
%%     ~~Describe \code{cluster.number} here~~
}
  \item{init_method}{
%%     ~~Describe \code{init_method} here~~
}
  \item{init_params}{
%%     ~~Describe \code{init_params} here~~
}
  \item{init_weights}{
%%     ~~Describe \code{init_weights} here~~
}
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
}
  \item{filtering}{
%%     ~~Describe \code{filtering} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
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
function (data, params, model = NULL, iter = 100, cluster.number = 4, 
    init_method = c("cluster", "random", "by"), init_params = NULL, 
    init_weights = NULL, mutation_position = "middle", filtering = c("single", 
        "unsolved"), trace = TRUE) 
{
    if (filtering == "single") {
        fi <- filter_unsolved_and_single(data = data, params = params)
        fdata <- fi$fdata
        fparams <- fi$fparams
    }
    else if (filtering == "unsolved") {
        fi <- filter_unsolved(data = data, params = params)
        fdata <- fi$fdata
        fparams <- fi$fparams
    }
    else {
        stop("Invalid argument filtering: must be either single or unsolved")
    }
    if (init_method == "cluster") {
        init <- initialize_by_clustering(data = fdata, params = fparams, 
            model = model, mutation_position = mutation_position, 
            cluster.number = cluster.number)
        iparameters <- init$iparameters
        iweights <- init$iweights
    }
    else if (init_method == "random") {
        init <- initialize_random(params = fparams, model = model, 
            cluster.number = cluster.number)
        iparameters <- init$iparameters
        iweights <- init$iweights
    }
    else if (init_method == "by") {
        if (is.null(init_params) || is.null(init_weights)) {
            stop("Chosen inititalization method requires init_params and init_weights arguments")
        }
        init <- initialize_by(init_params, init_weights, model = model, 
            cluster.number = cluster.number)
        iparameters <- init$iparameters
        iweights <- init$iweights
    }
    else {
        stop(paste("Invalid initialization method ", method))
    }
    print("Initial parameters:")
    print(iparameters)
    print(iweights)
    em_results <- em(data = fdata, parameters = iparameters, 
        model = model, weights = iweights, iter = iter, mutation_position = mutation_position, 
        cluster.number = cluster.number, trace = trace)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
