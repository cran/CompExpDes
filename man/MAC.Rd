\name{MAC}
\alias{MAC}
\title{
Maximum Absolute Correlation
}
\description{
Maximum Absolute Correlation (MAC) is the maximum absolute value rather than 1 of a correlation matrix.
}
\usage{
MAC(matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{
Input a matrix
}
}
\value{
It returns a maximum absolute value of the correlation matrix for a given matrix.
}

\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
library(CompExpDes)
lhd<-matrix(c(1,5,7,3,4,2,6,2,1,4,5,3,7,6,4,5,6,1,2,3,7),nrow=7,ncol=3,byrow=FALSE)
MAC(lhd)
}
\keyword{
maximum absolute correlation
}
\keyword{LHD}
\keyword{CompExpDes}
\keyword{Latin Hypercube Designs}
\keyword{Computer Experiments}
\keyword{Space-Filling Designs}