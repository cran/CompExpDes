\name{SLHDs}
\alias{SLHDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sliced Latin Hypercube Designs with Equal Size of Slices
}
\description{
This sliced LHDs are generated using a different new algorithm for a flexible parameter range which possesses a good space-filling property as whole design as well as for the slices.
}
\usage{
SLHDs(slices,factors,levels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slices}{
Number of slices(>=2)
}
\item{factors}{
Number of factors(>=2)
}
\item{levels}{
Number of levels(>=3)
}
}
\value{
Generates a Sliced LHD along with its parameters.
}
\author{
Ashutosh Dalal,
Cini varghese,
Rajender Parsad and Mohd Harun
}
\examples{
library(CompExpDes)
SLHDs(3,3,3)
}
\references{
Qian, P.Z.G. and Wu, C.F.J. (2009). Sliced space-filling designs. Biometrika, 96(4), 945–956. 
}
\keyword{LHD}
\keyword{CompExpDes}
\keyword{Latin Hypercube Designs}
\keyword{Computer Experiments}
\keyword{Space-Filling Designs}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
