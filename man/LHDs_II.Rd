\name{LHDs_II}
\alias{LHDs_II}
\title{
Latin Hypercube Designs (LHDs) for Any Numbers of Factors
}
\description{
For any number of factors, F (>2), this method will generate LHDs with runs ranges from F+2 to nC2, where n=2*F+1. Maxpro criterion measure, Phi_p measure also provided as a measure of space-filling and orthogonality measure maximum absolute correlation (MAC) value also provided.
}
\usage{
LHDs_II(levels,factors,weight,iterations)
}
\arguments{
  \item{levels}{ Ranges from (factors+2) to nC2, where n=2*factors+1}
  \item{factors}{ Any integer (>2)}
  \item{weight}{ Weight should be given to Maxpro, Phi_p and MAC such that sum is 1. Default it is 0.3, 0.3 and 0.4}
  \item{iterations}{ Number of iterations. By default it is 400.}
}
\value{
This function will provide a series of LHDs along with space-filling and orthogonality measures for any number.
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
\dontrun{
library(CompExpDes)
LHDs_II(20,3,c(0.4,0.2,0.4))
}
}
\references{
McKay, M.D., Beckman, R.J. and Conover, W.J. (1979). Comparison of three methods for selecting values of input variables in the analysis of output from a computer code. Technometrics, 21(2), 239-245.
}
\keyword{LHD}
\keyword{CompExpDes}
\keyword{Latin Hypercube Designs}
\keyword{Computer Experiments}
\keyword{Space-Filling Designs}