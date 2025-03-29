\name{wtLHDs_prime}
\alias{wtLHDs_prime}
\title{
Weighted Criteria-Based Latin Hypercube Designs (LHDs) for Prime Numbers
}
\description{
For prime number of factors, F (>2), this method will generate LHDs with levels, L ranges from F to F^2. Maxpro criterion measure, Phi_p measure also provided as a measure of space-filling and also as an orthogonality measure maximum absolute correlation (MAC) value also provided.
}
\usage{
wtLHDs_prime(levels,factors,weights,iterations)
}
\arguments{
  \item{levels}{ L,ranges from F to F^2}
  \item{factors}{ A prime number, F (>2)}
  \item{weights}{ Weight should be given to Maxpro, Phi_p and MAC such that sum is 1. Default it is 0.2, 0.4 and 0.4}
  \item{iterations}{ Number of iterations. By default it is 300.}
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
wtLHDs_prime(9,3,c(0.6,0,0.4))
}
\references{
McKay, M.D., Beckman, R.J. and Conover, W.J. (1979). Comparison of three methods for selecting values of input variables in the analysis of output from a computer code. Technometrics, 21(2), 239-245.
}
\keyword{LHD}
\keyword{CompExpDes}
\keyword{Latin Hypercube Designs}
\keyword{Computer Experiments}
\keyword{Space-Filling Designs}