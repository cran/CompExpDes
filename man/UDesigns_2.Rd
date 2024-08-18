\name{UDesigns_2}
\alias{UDesigns_2}
\title{
Uniform Designs with multiple factors
}
\description{
For v = n(n-1)/2, where n (>=5) is any odd number.  These are uniform designs interms of discrete discrepancy. It also provides number of factors, number of levels and number of runs and discrete discrepancy measure with its lower bound value.
}
\usage{
UDesigns_2(n)
}
\arguments{
  \item{n}{ any odd integer >=5
}
}
\value{
Returns a series of high dimensional uniform designs along with number of factors, levels, runs and discrete discrepancy measure with its lower bound value.
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
library(CompExpDes)
UDesigns_2(5)
}
\references{
Fang, K.T. (1980). The uniform design: application of number-theoretic methods in experimental design. Acta Math Appl Sin, 3, 363-372.
}
\keyword{Uniform designs}
\keyword{UD}
\keyword{Discrete Discrepancy}
\keyword{Multi dimensional uniform designs}
\keyword{LHD}
\keyword{CompExpDes}
\keyword{Latin Hypercube Designs}
\keyword{Computer Experiments}
\keyword{Space-Filling Designs}