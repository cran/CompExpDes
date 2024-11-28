\name{OLHDs_2F}
\alias{OLHDs_2F}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two Factor Orthogonal Latin Hypercube Designs
}
\description{
This OLHDs are generated using a new unique algorithm  which possesses a good space-filling property.
}
\usage{
OLHDs_2F(levels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{levels}{
Number of levels, 4t-1, where t =3,4,...
}
}
\value{
Generates two factor OLHDs along with its parameters.
}
\author{
Ashutosh Dalal,
Cini varghese,
Rajender Parsad and Mohd Harun
}
\examples{
library(CompExpDes)
OLHDs_2F(11)
}
\references{
Ye, K.Q.: Orthogonal column Latin hypercubes and their application in computer experiments. Journal of the American Statistical Association, 93(444), 1430–1439 (1998).  
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
