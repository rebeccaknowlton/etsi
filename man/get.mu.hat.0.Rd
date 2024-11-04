\name{get.mu.hat.0}
\alias{get.mu.hat.0}
\title{
Calculates the conditional expectation of the primary outcome from the control group in Study A, given the surrogate outcome, using kernel smoothing.
}
\description{
Calculates the conditional expectation of the primary outcome from the control group in Study A, given the surrogate outcome, using kernel smoothing.
}
\usage{
get.mu.hat.0(s, A.s0, A.y0)
}
\arguments{
  \item{s}{
s, the surrogate outcome
}
  \item{A.s0}{
A.s0, the vector of surrogate outcomes from the control group in Study A
}
  \item{A.y0}{
A.y0, the vector of primary outcomes from the control group in Study A
}
}
\value{
Returns the estimated conditional expectation of the primary outcome
}
\keyword{internal}