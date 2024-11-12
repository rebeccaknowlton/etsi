\name{etsi.design}
\alias{etsi.design}
\title{
Calculates the estimated power for a future Study B or the required sample size to achieve a desired level of power.
}
\description{
Function to help design a future Study B given data from Study A. Depending on the inputs, either calculates the estimated power in Study B, or calculates the required sample size in each treatment groupto achieve a desired level of power in Study B.
}
\usage{
etsi.design(Study.A, n.b0 = NULL, n.b1 = NULL, psi = NULL, w.range = NULL, 
            kappa = NULL, desired.power = NULL, iterations = 100)
}

\arguments{
  \item{Study.A}{
Study.A, a dataframe representing Study A containing the required columns A (treatment assignment), Y (primary outcome), S (surrogate outcome), and W (baseline covariate of interest). Optionally contains a column delta, which is an indicator for strong surrogacy.
}
  \item{n.b0}{
n.b0, sample size in Study B control group
}
  \item{n.b1}{
n.b1, sample size in Study B treatment group
}
  \item{psi}{
psi, the overall treatment effect you want to detect in Study B
}
  \item{w.range}{
w.range, a vector containing the lower and upper limits of the covariate W you plan to have in Study B
}
  \item{kappa}{
kappa, the threshold to determine strong surrogacy when the delta columns is not provided in Study.A. Surrogacy is considered sufficiently strong when the estimated PTE with respect to W is greater than kappa.
}
  \item{desired.power}{
The desired power in Study B, when you want to calculate required sample size
}
  \item{iterations}{
The number of iterations used for generalized cross-validation.
}
}
\value{
A list is returned:
  \item{power}{The estimated power in Study B; only given if desired.power = NULL.}
  \item{sample.size}{The required sample size in each treatment group in Study B; only given if desired.power is not NULL.}
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
  etsi.design(Study.A = exampledataA,
            n.b0 = 400, 
            n.b1 = 500, 
            iterations = 10)
}
