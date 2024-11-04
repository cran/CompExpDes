\name{Maxpro_Measure}
\alias{Maxpro_Measure}
\title{
Measure of Maxpro criterion
}
\description{
This function calculates Maxpro criterion for a given space-filling design. Lesser the value of it better the design, in the sense that the design has maximum spread in higher dimensional spaces.}
\usage{Maxpro_Measure(Design)}
\arguments{
  \item{Design}{
Provide a design in a matrix format
}
}
\value{
Provides Maxpro criterion value given by Joseph et al. (2015).
}
\references{
Joseph, V.R., Gul, E. and Ba, S. (2015). Maximum projection designs for computer experiments. Biometrika, 102 (2), 371-380.
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
Maxpro_Measure(lhd)
}
\keyword{Space-filling measure}
\keyword{maxpro}
\keyword{Computer experiments}
\keyword{LHD}
\keyword{CompExpDes}
\keyword{Latin Hypercube Designs}
\keyword{Space-Filling Designs}