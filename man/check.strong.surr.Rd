\name{check.strong.surr}
\alias{check.strong.surr}
\title{
Calculates whether each individual has a strong surrogate, based on PTE results from Study A and a desired threshold.
}
\description{
Calculates whether each individual has a strong surrogate, based on PTE results from Study A and a desired threshold.
}
\usage{
check.strong.surr(df, k, PTE.results)
}
\arguments{
  \item{df}{
The dataframe to determine strong surrogacy.
}
  \item{k}{
k, the threshold at which PTE > k qualifies as a strong surrogate
}
  \item{PTE.results}{
The PTE results (from the hetsurr::hetsurr.fun function) that should be used to determine strong surrogacy.
}
}
\value{
Returns a vector the same length as df, indicating for each row of df whether the surrogate is strong (1) or weak (0).
}
\keyword{internal}