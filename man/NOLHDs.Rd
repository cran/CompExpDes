\name{NOLHDs}
\alias{NOLHDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nearly Orthogonal Latin Hypercube Designs for Flexible Levels and Factors
}
\description{
This NOLHDs are generated using a new algorithm for a flexible parameter range which possesses a good space-filling property.
}
\usage{
NOLHDs(factors,levels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{factors}{
Number of factors(>=2)
}
\item{levels}{
Number of levels(>=factors+3)
}
}
\value{
Generates NOLHDs along with its parameters and maximum absolute correlation (MAC) value.
}
\author{
Ashutosh Dalal,
Cini varghese,
Rajender Parsad and Mohd Harun
}
\examples{
library(CompExpDes)
NOLHDs(2,5)
}
\references{
Cioppa, T.M. and Lucas, T.W.: Efficient nearly orthogonal and space-filling latin hypercubes. Technometrics, 49(1), 45–55 (2007).
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
