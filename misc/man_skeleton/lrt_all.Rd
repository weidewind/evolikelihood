\name{lrt_all}
\alias{lrt_all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lrt_all(mutation_position = "middle", fishy = TRUE, pack = "rootsolve", tag = "defaulttag", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
}
  \item{fishy}{
%%     ~~Describe \code{fishy} here~~
}
  \item{pack}{
%%     ~~Describe \code{pack} here~~
}
  \item{tag}{
%%     ~~Describe \code{tag} here~~
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
function (mutation_position = "middle", fishy = TRUE, pack = "rootsolve", 
    tag = "defaulttag", verbose = TRUE) 
{
    prots <- c("h1", "h3", "n1", "n2")
    for (pr in prots) {
        prot <- pr
        prot_data <- read.csv(paste(c(getwd(), "/input/", prot, 
            "_for_LRT.csv"), collapse = ""), stringsAsFactors = FALSE)
        splitted <- split(prot_data, list(prot_data$site, prot_data$ancestor_node), 
            drop = TRUE)
        if (verbose) {
            sink(file = paste(c(getwd(), "/output/", prot, "_all_parms_", 
                tag), collapse = ""))
        }
        prms <- parameters(splitted, mutation_position = mutation_position, 
            verbose = verbose, pack = pack, filter = FALSE)
        if (verbose) {
            sink()
        }
        lrt_procedure(data = splitted, prot = prot, tag = tag, 
            fishy = fishy, params = prms)
        sink()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
