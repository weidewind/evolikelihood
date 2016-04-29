\name{find_single_root}
\alias{find_single_root}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find_single_root(data, mutation_position, rkvector, jack = FALSE, pack = "rootsolve", verbose = TRUE)
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
function (data, mutation_position, rkvector, jack = FALSE, pack = "rootsolve", 
    verbose = TRUE) 
{
    if (pack != "rootsolve" && pack != "nleqslv") {
        stop("Incorrect argument pack: must be either 'rootsolve' or 'nleqslv'")
    }
    if (length(data) == 1) {
        mode <- "single"
        node_data <- data[[1]]
        rkvector <- c(1)
        names(rkvector) <- names(data)
        if (verbose) {
            print("Running in single mode")
        }
    }
    else {
        mode <- "group"
        if (verbose) {
            print("Running in group mode")
        }
        if (is.null(rkvector)) {
            stop("In group mode argument rkvector must be provided")
        }
    }
    if (mode == "group" || (mode == "single" && sum(node_data$event) != 
        0)) {
        if (verbose && mode == "single") {
            print(node_data)
            print(paste("site ", node_data[2, 1], " node ", node_data[2, 
                2]))
            print(paste("number of mutations ", sum(node_data$event)))
        }
        pars <- list(data = data, rkvector = rkvector, mutation_position = mutation_position)
        lambda_exp_root <- lambda_derivative_exp(pars)
        if (verbose) {
            print(paste("expon_lambda root ", lambda_exp_root))
        }
        p = -1
        init = 0.2
        count = 1
        while ((!is.na(p) && (p < 0 || (!is.null(solution_p$termcd) && 
            solution_p$termcd > 2) || (!is.null(solution_p$estim.precis) && 
            is.nan(solution_p$estim.precis)))) && count < 20) {
            print(paste(c("staring with p = ", init)))
            if (pack == "rootsolve") {
                if (jack) {
                  solution_p <- multiroot(f = p_derivative, start = c(init), 
                    jacfunc = p_derivative_jacfunc, jactype = "fullusr", 
                    parms = pars)
                }
                else {
                  solution_p <- multiroot(f = p_derivative, start = c(init), 
                    parms = pars)
                }
                p <- solution_p$root
            }
            else if (pack == "nleqslv") {
                solution_p <- nleqslv(fn = p_derivative, x = c(init), 
                  parms = pars)
                p <- solution_p$x
            }
            init <- init + 1
            count <- count + 1
        }
        if (pack == "rootsolve") {
            precision <- solution_p$estim.precis
            p_root <- solution_p$root
        }
        else {
            precision <- solution_p$fvec
            p_root <- solution_p$x
        }
        if (verbose) {
            print(paste("precision_p ", precision))
        }
        if (is.na(p_root)) {
            print(paste("no p_roots found"))
            c(p_root = NA, p_precision = NA, lambda_weib_root = NA, 
                lambda_exp_root = lambda_exp_root)
        }
        else if (is.na(precision)) {
            if (verbose) {
                print(paste("p_root is ", p_root, " but estimated precision for p is Na, won't try to estimate lambda"))
            }
            c(p_root = p_root, p_precision = NA, lambda_weib_root = NA, 
                lambda_exp_root = lambda_exp_root)
        }
        else {
            pars <- list(p = p_root, data = data, rkvector = rkvector, 
                mutation_position = mutation_position)
            lambda_weib_root <- lambda_derivative_weib(pars)
            if (verbose) {
                print(c(p_root = p_root, lambda_weib_root = lambda_weib_root))
            }
            c(p_root = p_root, p_precision = precision, lambda_weib_root = lambda_weib_root, 
                lambda_exp_root = lambda_exp_root)
        }
    }
    else {
        if (verbose && mode == "single") {
            print(paste("site ", node_data[2, 1], " node ", node_data[2, 
                2]))
            print("No mutations in the subtree, all roots NA")
        }
        c(p_root = NA, p_precision = NA, lambda_weib_root = NA, 
            lambda_exp_root = NA)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
