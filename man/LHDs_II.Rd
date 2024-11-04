\name{LHDs_II}
\alias{LHDs_II}
\title{
Latin Hypercube Designs (LHDs) for Any Numbers of Factors
}
\description{
For any even number of factors, F (>=4), this method will generate LHDs with levels, L ranges from F+2 to sC2, where s=2*F+1. Maxpro criterion measure, Phi_p measure also provided as a measure of space-filling and as an orthogonality measure, maximum absolute correlation (MAC) value also provided.
}
\usage{
LHDs_II(levels,factors,weight,iterations)
}
\arguments{
  \item{levels}{ Ranges from (F+2) to sC2, where s=2*F+1}
  \item{factors}{ F, any even number (>=4)}
  \item{weight}{ Weight should be given to Maxpro, Phi_p and MAC such that sum is 1. Default it is 0.3, 0.3 and 0.4}
  \item{iterations}{ Number of iterations. By default it is 400.}
}
\value{
This function will provide a series of LHDs along with space-filling and orthogonality measures for the generated LHDs.
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