\name{etsi.main}
\alias{etsi.main}
\title{
Estimates the pooled treatment effect quantity, standard error, and corresponding p-value.
}
\description{
Estimates the pooled treatment effect quantity in Study B, using the surrogate information learned from Study A. Provides the corresponding estimated standard error and p-value for testing the null hypothesis of a treatment effect equal to 0.
}
\usage{
etsi.main(Study.A, Study.B, kappa = NULL)
}

\arguments{
  \item{Study.A}{
Study.A, a dataframe representing Study A containing the required columns A (treatment assignment), Y (primary outcome), S (surrogate outcome), and W (baseline covariate of interest). Optionally contains a column delta, which is an indicator for strong surrogacy.
}
  \item{Study.B}{
Study.B, a dataframe representing Study B containing the required columns A (treatment assignment), Y (primary outcome), S (surrogate outcome), and W (baseline covariate of interest). Optionally contains a column delta, which is an indicator for strong surrogacy.
}
  \item{kappa}{
kappa, the threshold to determine strong surrogacy when the delta columns are not provided in Study.A and Study.B. Surrogacy is considered sufficiently strong when the estimated PTE with respect to W is greater than kappa.
}
}

\value{
A list is returned: 
  \item{delta.P}{The pooled treatment effect quantity for Study B}
  \item{se.delta.P}{Standard error of delta.P}
  \item{p.value}{p-value for testing the null hypothesis that delta.P = 0}
}
\references{
Knowlton, R., Parast, L. (2024) "Efficient Testing Using Surrogate Information." Under Review.
}
\author{
Rebecca Knowlton
}
\examples{
  data(exampledataA)
  data(exampledataB)
  names(exampledataA)
  names(exampledataB)
  etsi.main(Study.A = exampledataA,
            Study.B = exampledataB)
}