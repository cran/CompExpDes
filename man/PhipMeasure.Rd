\name{PhipMeasure}
\alias{PhipMeasure}
\title{
Phi_p criterion
}
\description{
For a given design Phi_p criterion is calculated. It is based on Morris and Mitchell (1995). When the designs are not in [0,1]^d form, lesser the value of Phi_p criterion better it is.
}
\usage{
PhipMeasure(design,p=15,q=2)
}
\arguments{
  \item{design}{
A design matrix is needed
}
  \item{p}{
Any positive integer. Default value of p = 15.
}
  \item{q}{
Any positive integer. Default value of q = 2. This implies that we are considering here Euclidean distance.
}
}
\value{
Generates Phi_p criterion value
}
\references{
Morris, M.D. and Mitchell, T.J. (1995). Exploratory designs for computer experiments. Journal of Statistical Planning and Inference, 43, 38-402. 
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and
Mohd Harun
}
\examples{
library(CompExpDes)
lhd<-matrix(c(1,5,7,3,4,2,6,2,1,4,5,3,7,6,4,5,6,1,2,3,7),nrow=7,ncol=3,byrow=FALSE)
PhipMeasure(lhd,p=15,q=2)
}
\keyword{Space-filling measure}
\keyword{phip}
\keyword{Computer experiments}

