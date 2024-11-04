\name{extrapolate.na}
\alias{extrapolate.na}
\title{
Function for extrapolating NA values from kernel smoothing
}
\description{
Function for extrapolating NA values from kernel smoothing
}
\usage{
extrapolate.na(y.strong, s.strong)
}
\arguments{
  \item{y.strong}{
y.strong, the vector of predicted primary outcome values for s.strong
}
  \item{s.strong}{
s.strong, the vector of strong surrogate information
}
}
\value{
Returns a vector y.strong where the NA values have been filled in with extrolated values based on s.strong
}
\keyword{internal}