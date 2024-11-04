\name{exampledataA}
\alias{exampledataA}
\docType{data}
\title{
Example Data (Study A)
}
\description{
Example Data (Study A)
}
\usage{data("exampledataA")}
\format{
  A data frame with 1000 observations from a treatment group and 1100 observations from a control group:
  \describe{
    \item{\code{A}}{the treatment assignment, where 1 indicates treatment and 0 indicates control}
    \item{\code{W}}{the baseline covariate of interest}
    \item{\code{S}}{the surrogate outcome}
    \item{\code{Y}}{the primary outcome}
    \item{\code{delta}}{the indicator for strong surrogacy, where 1 indicates strong and 0 indicates weak}
  }
}
\examples{
data(exampledataA)
names(exampledataA)
}
