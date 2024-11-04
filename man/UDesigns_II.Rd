\name{UDesigns_II}
\alias{UDesigns_II}
\title{
Uniform Designs with Multiple Factors with Minimal Runs
}
\description{
 In this series, the Uniform Designs (UDs) are high dimensional with lesser number of runs will always attain lower bound of discrete discrepancy. They are available for any even number of factors, F (>= 4) with F(F+1) levels each.
}
\usage{
UDesigns_II(factors)
}
\arguments{
  \item{factors}{ any even number >=4
}
}
\value{
Returns a series of high dimensional UDs along with number of factors, levels, runs, MAC value and discrete discrepancy measure along with its lower bound value.
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
library(CompExpDes)
UDesigns_II(4)
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