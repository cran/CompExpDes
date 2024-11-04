\name{Best_Model}
\alias{Best_Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find Best Model 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function will try to find out a significant model for each combinations based on adjusted R^2. Then user need to select which model they want to use.
}
\usage{
Best_Model(model, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
Provide a vector that contains all the individual terms present in a full model
}
  \item{data}{
Provide data in a matrix or data frame format where you want to fit the model
}
}
\value{
Generate a list of significant models for various combinations of factors.
}
\author{
Ashutosh Dalal,
Cini Varghese,
Rajender Parsad and Mohd Harun
}
\examples{
\dontrun{
library(CompExpDes)
# Sample data
data <- data.frame(
  x1 = c(1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 5.4),
  x2 = c(50, 25, 5, 30, 55, 45, 20, 10, 35, 60, 40, 15),
  x3 = c(2.5, 6.0, 4.0, 1.0, 5.5, 4.5, 3.0, 2.0, 6.5, 3.5, 1.5, 5.0),
  x4 = c(45, 25, 55, 35, 65, 15, 70, 20, 50, 30, 60, 40),
  y = c(0.0795, 0.0118, 0.0109, 0.0991, 0.1266, 0.0717, 0.1319, 0.0900, 0.1739, 
  0.1176, 0.1836, 0.1424)
)
# List of terms in the polynomial model
model <- list('x1', 'x2', 'x3', 'x4', 'x1:x2', 'x1:x3', 'x1:x4', 
                   'x2:x3', 'x2:x4', 'x3:x4', 'I(x1^2)', 
                   'I(x2^2)', 'I(x3^2)', 'I(x4^2)')

Best_Model(model,data)
}
}

