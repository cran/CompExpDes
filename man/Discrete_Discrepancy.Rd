\name{Discrete_Discrepancy}
\alias{Discrete_Discrepancy}
\title{
Measure of Discrete Discrepancy
}
\description{
Discrete Discrepancy is a measure of uniformity for any uniform design. Lesser the value of Discrete Discrepancy measure, better is the uniform design. 
}
\usage{
Discrete_Discrepancy(Design,a,b)
}
\arguments{
  \item{Design}{
A matrix
}
  \item{a}{
Any value a>b>0. By default it is set to 1.
}
  \item{b}{
Any value a>b>0. By default it is set to 0.5.
}
}
\value{
The function calculates the value of Discrete Discrepancy measure and its lower bound for a given design.
}
\references{
Qin H, Fang KT (2004)<DOI:10.1007/s001840300296> Discrete discrepancy in factorial designs. Metrika, 60, 59-72.
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
library(CompExpDes)
lhd1<-matrix(c(1,5,7,3,4,2,6,2,1,4,5,3,7,6,4,5,6,1,2,3,7),nrow=7,ncol=3,byrow=FALSE)
lhd2<-cbind(lhd1[,3],lhd1[,2],lhd1[,1])
ud<-rbind(lhd1,lhd2)
Discrete_Discrepancy(ud, 1, 0.5)
}
\keyword{
Discrete Discrepancy
}
\keyword{
Uniformity measure
}
\keyword{LHD}
\keyword{CompExpDes}
\keyword{Latin Hypercube Designs}
\keyword{Computer Experiments}
\keyword{Space-Filling Designs}