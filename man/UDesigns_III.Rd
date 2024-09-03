\name{UDesigns_III}
\alias{UDesigns_III}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Orthogonal Uniform Designs for Two and Four Factors (Even number v)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
For even number v (>=6) this function will provide total three uniform designs of two factors (2) and four factors (1).
}
\usage{
UDesigns_III(v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
%%     ~~Describe \code{x} here~~
An even number of factors,v >=6
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Uniform Design 1 and Uniform Design 3 are orthogonal/ nearly orthogonal but Uniform Design 2 is always orthogonal design.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
This function will generate 3 Uniform Designs along with the number of levels,factors, runs and MAC value.
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
UDesigns_III(10)
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