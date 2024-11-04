\name{Meeting_Number}
\alias{Meeting_Number}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Coincidence (or Meeting) numbers between rows 
}
\description{
Finding out Maximum coincidence (or Meeting) number between unique pair of rows. 
}
\usage{
Meeting_Number(matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{
Provide any matrix
}
}
\value{
This function provides the maximum coincidence number between any pair of rows of for a given matrix.
}
\author{
  Ashutosh Dalal,
  Cini Varghese,
  Rajender Parsad and
  Mohd Harun
}
\examples{
library(CompExpDes)
  mat<-matrix(c(1,2,3,3,2,1,4,2,1),nrow=3,byrow=TRUE)
  Meeting_Number(mat)
}
\keyword{Meeting number}
\keyword{Coincidence number}
\keyword{CompExpDes}
