\name{MaxproLHD_1}
\alias{MaxproLHD_1}
\title{
Maxpro Latin Hypercube Designs (LHDs) for Prime Numbers
}
\description{
For a prime number s (>2), this method will generate an optimal maxpro LHD with runs v = s^2. Maxpro criterion measure, Phi_p measure also provided as a measure of space-filling and orthogonality measure maximum absolute correlation value also provided.

}
\usage{
MaxproLHD_1(prime_number)
}
\arguments{
  \item{prime_number}{ A prime number (>2)
}
}
\value{
This function will provide a series of optimal LHDs along with space-filling and orthogonality measures for a prime numbers.
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
library(CompExpDes)
MaxproLHD_1(5)
}
\references{
McKay, M.D., Beckman, R.J. and Conover, W.J. (1979). Comparison of three methods for selecting values of input variables in the analysis of output from a computer code. Technometrics, 21(2), 239-245.
}
\keyword{Maxpro LHD}
\keyword{LHD}
\keyword{Maxpro designs}
\keyword{Maximum projection designs}
\keyword{Latin Hypercube Designs}