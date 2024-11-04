\name{UDesigns_III}
\alias{UDesigns_III}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nearly Orthogonal Uniform Designs for Two and Four Factors 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function will provide nearly orthogonal uniform designs (UDs) for number of factors, F = 2 and 4 but a flexible number of levels, L >=3.
}
\usage{
UDesigns_III(levels,factors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{levels}{
%%     ~~Describe \code{x} here~~
Number of levels, L >=3
}
  \item{factors}{
%%     ~~Describe \code{x} here~~
Number of factors 2 or 4
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
This function will generate 3 Uniform Designs along with the number of levels,factors, runs, MAC value and discrete discrepancy value along with its lower bound value.
}
\references{
%% ~put references to the literature/web site here ~
Fang, K.T. (1980). The uniform design: application of number-theoretic methods in experimental design. Acta Math Appl Sin, 3, 363-372.
}
\author{
Ashutosh Dalal,
Cini varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
\dontrun{
library(CompExpDes)
UDesigns_III(3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
\keyword{Uniform designs}
\keyword{UD}
\keyword{Discrete Discrepancy}
\keyword{Orthogonal uniform designs}
\keyword{Two factor orthogonal uniform designs}
\keyword{CompExpDes}
\keyword{Computer Experiments}
\keyword{Space-Filling Designs}