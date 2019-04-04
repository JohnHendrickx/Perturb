\name{reclassify}
\alias{reclassify}
\alias{print.reclassify}
\title{Called by perturb to calculate reclassification tables}
\description{
  \code{reclassify} is called by \code{\link{perturb}} to calculate reclassification probabilities for categorical variables. Use separately to experiment with reclassification probabilities.
}
\usage{
reclassify(varname, pcnt = NULL, adjust = TRUE, bestmod = TRUE,
min.val = .1, diag = NULL, unif = NULL, dist = NULL, assoc = NULL, data=NULL)

\method{print}{reclassify}(x, dec.places = 3, full = FALSE, ...)
}
\arguments{
  \item{varname}{a factor to be reclassified}
  \item{pcnt}{initial reclassification percentages}
  \item{adjust}{makes the expected frequency distribution of the reclassified variable equal to that of the original}
  \item{bestmod}{imposes an appropriate pattern of association between the original and the reclassified variable}
  \item{min.val}{value to add to empty cells of the initial expected table when estimating the best model}
  \item{diag}{The odds of same versus different category reclassification}
  \item{unif}{Controls short distance versus long distance reclassification for ordered variables}
  \item{dist}{alternative parameter for short versus long distance reclassification}
  \item{assoc}{a matrix defining a loglinear pattern of association}
  \item{data}{the data frame containing \code{varname}}
  \item{x}{a \code{reclassify} object to be printed}
  \item{dec.places}{number of decimal places to use when printing}
  \item{full}{if TRUE, some extra information is printed}
  \item{\dots}{arguments to be passed on to or from other methods. Print options for class \code{matrix} may be used, e.g. \code{print.gap}}
}

\details{
\code{reclassify} creates a table of reclassification probabilities for \emph{varname}. By default, the reclassification probabilities are defined so that the expected frequency distribution of the reclassified variable is identical to that of the original. In addition, a meaningful pattern of association is imposed between the original and the reclassified variable.\code{reclassify} is called by \code{\link{perturb}} to calculate reclassification probabilities for categorical variables. \code{reclassify} can be used separately to find a suitable reclassification probabilities.

\code{Reclassify} has several options but the most relevant will generally be the \code{pcnt} option. The argument for \code{pcnt} can be

\itemize{
	\item a scalar
	\item a vector of length \eqn{n}
	\item a vector of length \eqn{n^2}, where \eqn{n} is the number of categories of the variable to be reclassified.
}

If the argument for \code{pcnt} is a scalar, its value is taken to be the percentage of cases to be reclassified to the same category, which is the same for all categories. A table of initial reclassification probabilities for the original by the reclassified variable is created with this value divided by 100 on the diagonal and equal values on off-diagonal cells.

If the argument for \code{pcnt} is a vector of length \eqn{n}, its values indicate the percentage to be reclassified to the same category for each category separately. These values divided by 100 form the diagonal of the table of initial reclassification probabilities. Off-diagonal cells have the same values for rows so that the row sum is equal to 1.

If the argument for \code{pcnt} is a vector of length \eqn{n^2}, its values form the table of initial reclassification probabilities. \code{\link{prop.table}} is used to ensure that these values sum to 1 over the columns. Specifying a complete table of initial reclassification probabilities will be primarily useful when an ordered variable is being reclassified.

\code{Reclassify} prints an initial table of reclassification probabilities based on the \code{pcnt} option. This table is not used directly though but \emph{adjusted} to make the expected frequencies of the reclassified variable identical to those of the original. In addition, a meaningful pattern of association is imposed between the original and the reclassified variable. Details are given in the section \emph{\dQuote{Adjusting the reclassification probabilities}}.

Knowledgeable users can specify a suitable pattern of association directly, bypassing the pcnt option. Details are given in the section
\emph{\dQuote{Specifying a pattern of association directly}}.
}

\value{
  An object of class \code{reclassify}. By default, \code{print.reclassify} prints the variable name and the \code{reclass.prob}. If the \code{full} option is used with \code{print.reclassify}, additional information such as the initial reclassification probabilities, initial expected table, best model, are printed as well.

  \item{variable}{The variable specified}
  \item{reclass.prob}{Row-wise proportions of \code{fitted.table}}
  \item{cum.reclass.prob}{Cumulative row-wise proportions}
  \item{exptab$init.pcnt}{initial reclassification probabilities  (option \code{pcnt})}
  \item{exptab$init.tbl}{initial expected frequencies (option \code{pcnt})}
  \item{bestmod}{The best model found for the table of initial expected frequencies (option \code{pcnt})}
  \item{assoc}{The log pattern of association specified using \code{pcnt} and \code{bestmod=FALSE}}
  \item{coef}{The coefficients of a fitted loglinear model}
  \item{fitted.table}{The adjusted table of expected frequencies}
}

\references{
Goodman, Leo A. (1984). \cite{The analysis of cross-classified data having ordered categories.} Cambridge, Mass.: Harvard University Press.

Hendrickx, J. (2004). Using standardised tables for interpreting loglinear models. \cite{Quality & Quantity} 38: 603-620.

Hendrickx, John, Ben Pelzer. (2004). \cite{Collinearity involving ordered and unordered categorical variables}. Paper presented at the RC33 conference in Amsterdam, August 17-20 2004. Available at \url{https://www.researchgate.net/publication/235994590_Collinearity_involving_ordered_and_unordered_categorical_variables}

Hout, M. (1983). \cite{Mobility tables}. Beverly Hills: Sage Publications.

Kaufman, R.L., & Schervish, P.G. (1986). Using adjusted crosstabulations to interpret log-linear relationships. \cite{American Sociological Review} 51:717-733
}
\author{John Hendrickx \email{John_Hendrickx@yahoo.com}}

\section{Adjusting the reclassification probabilities}{
A problem with the initial reclassification probabilities created using \code{pcnt} is that the expected frequencies of the reclassified variable will not be the same as those of the original. Smaller categories will become larger in the expected frequencies, larger categories will become smaller. This can be seen in the column marginal of the initial table of expected frequencies in the \code{reclassify} output. This could have a strong impact on the standard errors of reclassified variables, particularly as categories differ strongly in size.

To avoid this, the initial expected table is \emph{adjusted} so that the column margin is the same as the row margin, i.e. the expected frequencies of the reclassified variable are the same as those of the original. Use \code{adjust=FALSE} to skip this step. In that case the initial reclassification probabilities are also the final reclassification probabilities.

A second objection to the initial reclassification probabilities is that the pattern of association between the original and the reclassified variable is arbitrary. The association between some combinations of categories is higher than for others. \code{Reclassify} therefore derives an appropriate pattern of association for the initial expected table of the original by reclassified variable. This pattern of association is used when \dQuote{adjusting} the marginals to make the frequency distribution of the reclassified variable identical to that of the original. Use the option \code{bestmod=FALSE} to skip this step.

The patterns of association used by reclassify are drawn from loglinear models for square tables, also known as \dQuote{mobility models} (Goodman 1984, Hout 1983). Many texts on loglinear modelling contain a brief discussion of such models as well. For unordered variables, a \dQuote{quasi-independent} pattern of association would be appropriate. Under quasi-independent association, the row variable is independent of the column variable if the diagonal cells are ignored.

If the argument for \code{pcnt} was a scalar, \code{reclassify} fits a \dQuote{quasi-independent (constrained)} model. This model has a single parameter \code{diag} which indicates the log-odds of same versus different reclassification. This log-odds is the same for all categories. If the argument was of vector of length \eqn{n}, then a regular quasi-independence model is fitted with parameters \code{diag1} to \code{diag}\emph{n}. These parameters indicate the log-odds of same versus different category reclassification, which is different for each category. For both models, the reclassified category is independent of the original category if the diagonal cells are ignored.

If the argument for \code{pcnt} was a vector of length \eqn{n^2}, \code{reclassify} fits two models, a \dQuote{quasi-distance model} and a \dQuote{quasi-uniform association} model, and selects the one with the best fit to the initial expected table. Both have the \code{diag} parameter of the \dQuote{quasi-independence (constrained)} model. An additional parameter is added to make short distance reclassification more likely than long distance reclassification. The quasi-uniform model is stricter: it makes reclassification less likely proportionately to the squared difference between the two categories. The distance model makes reclassification less likely proportionately to the absolute difference between the two categories.

In some cases, the initial expected table based on the \code{pcnt} option contains empty cells. To avoid problems when estimating the best model for this table, a value of .1 is added to these cells. Use the \code{min.val} option to specify a different value.
}

\section{Specifying a pattern of association directly}{
If the \code{pcnt} option is used, \code{reclassify} automatically determines a suitable pattern of association between the original and the reclassified variable. Knowledgeable users can also specify a pattern of association directly. The final reclassification probabilities will then be based on these values. Built-in options for specifying the loglinear parameters of selected mobility models are:

\describe{
	\item{diag}{quasi-independence constrained  (same versus different category reclassification)}
	\item{unif}{uniform association (long versus short distance reclassification for ordered categories)}
	\item{dist}{linear distance model (allows more long distance reclassification than uniform association)}
}

The \code{assoc} option can be used to specify an association pattern of one's own choice. The elements of \code{assoc} should refer to matrices with an appropriate loglinear pattern of association. Such matrices can be created in many ways. An efficient method is:

\code{wrk<-diag(table(}\emph{factor}\code{))}\cr
\code{myassoc<-abs(row(wrk)-col(wrk))*-log(5)}

This creates a square diagonal matrix called \code{wrk} with the same number of rows and columns as the levels of \emph{factor}. \code{row(wrk)} and \code{col(wrk)} can now be used to define a loglinear pattern of association, in this case a distance model with parameter 5. \code{reclassify} checks the length of the matrix equals \eqn{n^2}, where \emph{n} is the number of categories of \code{varname} and ensures that the pattern of association is symmetric.
}

\section{Imposing a pattern of association}{
A table with given margins and a given pattern of association can be created by
\itemize{
\item estimating a loglinear model of independence for a table with the desired margins
\item while specifying the log pattern of association as an offset variable (cf. Kaufman & Schervish (1986), Hendrickx (2004).
}
The body of the table is unimportant as long as it has the appropriate margins. The predicted values of the model form a table with the desired properties.

The expected table of the original by the reclassified variable is adjusted by creating a table with the frequency distribution of the original variable on the diagonal cells. This table then has the same marginals for the row and column variables. The pattern of association is determined by the reclassify options. If \code{pcnt}  is used and \code{bestmod=TRUE} then the predicted values of the best model are used as the offset variable. If \code{bestmod=FALSE}, the log values of the initial expected table are made symmetric and used as the offset variable. If a loglinear model was specified directly, a variable is created in the manner of the \code{assoc} example.

A small modification in procedure is that reclassify uses a model of equal main effects rather than independence. Since the pattern of association is always symmetric, the created table will then also be exactly symmetric with the frequency distribution of the original variable as row and column marginal.
}

\seealso{
\code{\link{perturb}}, \code{\link{colldiag}}, \code{[car]}\code{\link[car]{vif}}, \code{[rms]}\code{\link[rms]{vif}}
}
\examples{
library(car)
data(Duncan)

reclassify("type",pcnt=95,data=Duncan)
}
\keyword{regression}
