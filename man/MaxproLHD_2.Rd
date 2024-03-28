\name{MaxproLHD_2}
\alias{MaxproLHD_2}
\title{
Maxpro Latin Hypercube Designs (LHDs) for odd numbers
}
\description{
For an odd number s, this method will generate 3 Maxpro-LHDs with runs v = s(s-1)/2.
Maxpro criterion measure, Phi_p measure also provided as a measure of space-filling and orthogonality measure maximum absolute correlation value also provided.
}
\usage{
MaxproLHD_2(s)
}
\arguments{
  \item{s}{ An odd number (s>=5)
}
}
\value{
This function will provide a series of 3 Maxpro-LHDs along with space-filling and orthogonality measures for an odd numbers.
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
library(CompExpDes)
MaxproLHD_2(7)
}
\references{
McKay, M.D., Beckman, R.J. and Conover, W.J. (1979). Comparison of three methods for selecting values of input variables in the analysis of output from a computer code. Technometrics, 21(2), 239-245.
}
\keyword{Maxpro LHD}
\keyword{LHD}
\keyword{Maxpro designs}
\keyword{Maximum projection designs}
\keyword{Latin Hypercube Designs}