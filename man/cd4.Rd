\name{cd4}
\docType{data}
\alias{cd4}
\title{Sample Data from "Applied Longitudinal Analysis"}
\description{
  Sample data from Chapter 8 of "Applied Longitudinal Analysis" by Garrett Fitzmaurice, Nan Laird & James Ware.
  The data are from the "ACTG Study 193A". This was a randomised, double-blind study of 1309 AIDS patients where subjects with dual combinations of HIV-1 inhibitors were compared to subjects with triple combinations of HIV-1 reverse transciptase inhibitors. The dependent variable in the dataset is logcd4, calculated as the log of the CD4 counts plus 1. The CD4 counts were assessed at baseline and between 1 and 9 times during the 40 weeks follow-up. The timing of assessments relative to baseline is measured in weeks.
}

\usage{data(cd4)}
\format{
  A data frame containing the following 6 variables.
  \describe{
    \item{id}{Subject id}
    \item{group}{Treatment Group}
    \item{age}{Age at baseline (years)}
    \item{sex}{Sex (1=male, 0=female)}
    \item{trt}{Treatment (1=Triple therapy, 0=Dual therapy)}
    \item{week}{Time since baseline (weeks)}
    \item{week_16}{Time since week 16 (max of 0, week - 16)}
    \item{w16}{week_16 - week}
    \item{logcd4}{log CD4 count plus 1}
  }
}
\source{https://content.sph.harvard.edu/fitzmaur/ala2e/}
\references{
Garrett Fitzmaurice, Nan Laird & James Ware. (2004).
\cite{Applied Longitudinal Analysis}. New York: John Wiley & Sons.
}
\examples{
data(cd4)
}
\keyword{datasets}
