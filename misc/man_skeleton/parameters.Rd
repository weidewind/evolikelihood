\name{parameters}
\alias{parameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parameters(data, mutation_position = "middle", rkvector, filter = TRUE, jack = FALSE, pack = "rootsolve", verbose = FALSE)
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
  \item{filter}{
%%     ~~Describe \code{filter} here~~
}
  \item{jack}{
%%     ~~Describe \code{jack} here~~
}
  \item{pack}{
%%     ~~Describe \code{pack} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (data, mutation_position = "middle", rkvector, filter = TRUE, 
    jack = FALSE, pack = "rootsolve", verbose = FALSE) 
{
    ps <- lapply(names(data), function(elm, mut_pos) {
        mutation_position <- mut_pos
        node_data <- data[elm]
        node_roots <- find_single_root(node_data, mutation_position, 
            rkvector, jack = jack, pack = pack, verbose = verbose)
        if (!filter || filter && !is.na(node_roots) && all(is.finite(node_roots)) && 
            node_roots["p_precision"] < 1e-05) {
            c(node = elm, lambda_exp = node_roots["lambda_exp_root"], 
                lambda_weib = node_roots["lambda_weib_root"], 
                p = node_roots["p_root"], p_precision = node_roots["p_precision"], 
                events = sum(node_data[[1]]$event))
        }
    }, mut_pos = mutation_position)
    if (filter) {
        ps <- Filter(Negate(is.null), ps)
    }
    params <- data.frame(matrix(unlist(ps), nrow = length(ps), 
        byrow = T), stringsAsFactors = FALSE)
    names(params) <- c("node", "lambda_exp_root", "lambda_weib_root", 
        "p_root", "p_precision", "events")
    params <- transform(params, lambda_exp_root = as.numeric(lambda_exp_root), 
        lambda_weib_root = as.numeric(lambda_weib_root), p_root = as.numeric(p_root), 
        p_precision = as.numeric(p_precision), events = as.numeric(events))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
