\name{supplier}
\alias{supplier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
supplier(data, rkvector, elm, mutation_position)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{rkvector}{
%%     ~~Describe \code{rkvector} here~~
}
  \item{elm}{
%%     ~~Describe \code{elm} here~~
}
  \item{mutation_position}{
%%     ~~Describe \code{mutation_position} here~~
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
function (data, rkvector, elm, mutation_position) 
{
    node_data <- data[[elm]]
    r <- rkvector[elm]
    if (mutation_position != "start" && mutation_position != 
        "end" && mutation_position != "middle") {
        stop("Value of mutation_position parameter must be either \"start\" or \"end\" or \"middle\"")
    }
    if (mutation_position == "start") {
        no_events_branches = node_data[node_data["event_indicator"] == 
            0, ]
    }
    else if (mutation_position == "end") {
        no_events_branches = node_data
    }
    else if (mutation_position == "middle") {
        no_events_branches <- node_data[node_data["event_indicator"] == 
            0, ]
        first_halves = node_data[node_data["event_indicator"] == 
            1, ]
        first_halves[, "t_branch_end"] <- (first_halves[, "t_branch_end"] + 
            first_halves[, "t_branch_start"])/2
        no_events_branches = rbind(no_events_branches, first_halves)
    }
    if (mutation_position == "middle") {
        second_halves = node_data[node_data["event_indicator"] == 
            1, ]
        second_halves[, "t_branch_start"] <- (second_halves[, 
            "t_branch_end"] + second_halves[, "t_branch_start"])/2
        events_branches = second_halves
    }
    else {
        events_branches = node_data[node_data["event_indicator"] == 
            1, ]
    }
    list(r = r, no_events_branches = no_events_branches, events_branches = events_branches)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
