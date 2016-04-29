\name{lrt_procedure}
\alias{lrt_procedure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lrt_procedure(data, prot, tag, fishy = FALSE, threshold = 3.84, mutation_position = "end", pack = "rootsolve", params = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{prot}{
%%     ~~Describe \code{prot} here~~
}
  \item{tag}{
%%     ~~Describe \code{tag} here~~
}
  \item{fishy}{
%%     ~~Describe \code{fishy} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
}
  \item{pack}{
%%     ~~Describe \code{pack} here~~
}
  \item{params}{
%%     ~~Describe \code{params} here~~
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
function (data, prot, tag, fishy = FALSE, threshold = 3.84, mutation_position = "end", 
    pack = "rootsolve", params = NULL) 
{
    sink(file = paste(c(getwd(), "/output/", prot, "_single_root_", 
        tag), collapse = ""), append = FALSE, type = c("output", 
        "message"), split = FALSE)
    lratios <- lapply(names(data), function(elm, mutation_position) {
        mutation_position <- mutation_position
        node_data <- data[elm]
        if (is.null(params)) {
            node_roots <- as.list(find_single_root(node_data, 
                mutation_position, pack = pack))
        }
        else {
            node_roots <- params[params$node == elm, ]
        }
        if (!is.na(node_roots) && all(!is.na(node_roots)) && 
            node_roots$p_precision < 1e-05) {
            lr <- lrt(node_data, lambda_exp = node_roots$lambda_exp_root, 
                lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, 
                fishy = TRUE, verbose = FALSE)
            row <- c(node = elm, lambda_exp = node_roots$lambda_exp_root, 
                lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, 
                p_precision = node_roots$p_precision, lr = lr)
            if (!is.na(lr) && lr > threshold) {
                print(c(node = elm, lambda_exp = node_roots$lambda_exp_root, 
                  lambda_weib = node_roots$lambda_weib_root, 
                  p = node_roots$p_root, p_precision = node_roots$p_precision, 
                  lr = lr))
            }
        }
        else {
            row <- c(node = elm, lambda_exp = node_roots$lambda_exp_root, 
                lambda_weib = node_roots$lambda_weib_root, p = node_roots$p_root, 
                p_precision = node_roots$p_precision, lr = NA)
        }
    }, mutation_position)
    sink()
    lratios
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
