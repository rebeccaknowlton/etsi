\name{exampledataB}
\alias{exampledataB}
\docType{data}
\title{
Example Data (Study B)
}
\description{
Example Data (Study B)
}
\usage{data("exampledataB")}
\format{
  A data frame with 500 observations from a treatment group and 400 observations from a control group:
  \describe{
    \item{\code{A}}{the treatment assignment, where 1 indicates treatment and 0 indicates control}
    \item{\code{W}}{the baseline covariate of interest}
    \item{\code{S}}{the surrogate outcome}
    \item{\code{Y}}{the primary outcome}
    \item{\code{delta}}{the indicator for strong surrogacy, where 1 indicates strong and 0 indicates weak}
  }
}
\examples{
data(exampledataB)
names(exampledataB)
}