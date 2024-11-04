\name{calculate.se}
\alias{calculate.se}
\title{
Estimates the standard error of the pooled treatment effect quantity.
}
\description{
Estimates the standard error of the pooled treatment effect quantity.
}
\usage{
calculate.se(y1.weak, y0.weak, y1.strong, y0.strong)
}
\arguments{
  \item{y1.weak}{
y1.weak, the primary outcome in the treatment group with weak surrogacy
}
  \item{y0.weak}{
y0.weak, the primary outcome in the control group with weak surrogacy
}
  \item{y1.strong}{
y1.strong, the predicted primary outcome given the surrogate in the treatment group with strong surrogacy
}
  \item{y0.strong}{
y0.strong, the predicted primary outcome given the surrogate in the control group with strong surrogacy
}
}
\value{
Returns the estimated standard error.
}
\keyword{internal}