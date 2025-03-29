\name{wtLHDs}
\alias{wtLHDs}
\title{
Weighted Criteria-Based Latin Hypercube Designs (LHDs) for Any Numbers of Factors (>=2)
}
\description{
For any number of factors, F (>=2), this method will generate LHDs with levels, L ranges from F to sC2, where s=2*F+1. Maxpro criterion measure, Phi_p measure also provided as a measure of space-filling and as an orthogonality measure, maximum absolute correlation (MAC) value also provided.
}
\usage{
wtLHDs(levels,factors,weights,iterations,population)
}
\arguments{
  \item{levels}{ Ranges from F to sC2, where s=2*F+1}
  \item{factors}{ F, any number (>=2)}
  \item{weights}{ Weight should be given to MAC, Phi_p and Maxpro such that sum is 1. Default it is 0.2, 0.4 and 0.4}
  \item{iterations}{ Number of iterations. By default it is 300}
\item{population}{ Population size. By default it is 100}
}
\value{
This function will provide a series of weighted criteria-based LHDs along with space-filling and orthogonality measures for the generated LHDs.
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
library(CompExpDes)
wtLHDs(20,3,c(0.4,0.2,0.4))
}
\references{
McKay, M.D., Beckman, R.J. and Conover, W.J. (1979). Comparison of three methods for selecting values of input variables in the analysis of output from a computer code. Technometrics, 21(2), 239-245.
}
\keyword{LHD}
\keyword{CompExpDes}
\keyword{Latin Hypercube Designs}
\keyword{Computer Experiments}
\keyword{Space-Filling Designs}