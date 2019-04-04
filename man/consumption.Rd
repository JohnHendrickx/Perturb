\name{consumption}
\docType{data}
\alias{consumption}
\title{The Consumption Function}
\description{
  Example from pp 149-154 of Belsley (1991), Conditioning Diagnostics
}
\usage{data(consumption)}
\format{
  A data frame with 28 observations on the following 5 variables.
  \describe{
    \item{year}{1947 to 1974}
    \item{cons}{total consumption, 1958 dollars}
    \item{rate}{the interest rate (Moody's Aaa)}
    \item{dpi}{disposable income, 1958 dollars}
    \item{d_dpi}{annual change in disposable income}
  }
}
\source{pp 149-154 of Belsley (1991), Conditioning Diagnostics}
\references{
Belsley, D.A. (1991).
\cite{Conditioning diagnostics, collinearity and weak data in regression}. New York: John Wiley & Sons.
}
\examples{
data(consumption)
}
\keyword{datasets}
