\name{perturb}
\alias{perturb}
\alias{summary.perturb}
\alias{print.summary.perturb}
\alias{GetData.default.perturb }
\alias{GetData.lme.perturb }
\alias{GetNames.default.perturb}
\alias{GetNames.lme.perturb }
\alias{coefs.default.perturb }
\alias{coefs.multinom.perturb}
\alias{coefs.merMod.perturb }
\alias{coefs.lme.perturb }
\alias{Repeatedcoefs.default.perturb}
\alias{Repeatedcoefs.lme.perturb }
\alias{VarCorr2.default.perturb}
\alias{VarCorr2.merMod.perturb }
\alias{VarCorr2.lme.perturb}
\alias{GetWarnings.default.perturb }
\alias{GetWarnings.merMod.perturb}
\title{Perturbation analysis to evaluate collinearity}
\description{
 Adds random noise to selected variables to evaluate collinearity. Suitable for other models than linear regression and for categorical independent variables.
 Version 3 can add add random noise to random and fixed effects in mixed models estimated by \code{\link[lme4]{lmer}} and \code{\link[nlme]{lme}}
}
\usage{
perturb(mod, pvars = NULL, prange = NULL, ptrans = NULL, postrun=NULL,
pfac = NULL, uniform = FALSE, niter = 100)

\method{summary}{perturb}(object,dec.places=3,full=FALSE,...)

\method{print.summary}{perturb}(x,\dots)

}
\arguments{
  \item{mod}{
  A model object, not necessarily type lm
  }
  \item{pvars}{
   Contains an array of variables to be perturbed. Random values are added to the variable, after which the model is re-estimated.
   }
  \item{prange}{
   Contains an array of values determining the magnitude of perturbations. There should be as many \emph{prange} values as \emph{pvars} variables.
  }
  \item{ptrans}{
  Contains an array of transformations to be applied. Each element should contain valid \R syntax in quotes for deriving one of the \emph{pvars} as a function of other variables.
  }
  \item{pfac}{
  Contains a list of categorical variables to be perturbed and parameters controlling the reclassification process. The first component must be a factor name. The name for the first component is ignored. Other list components should correspond with options for \code{\link{reclassify}}. The usage of these parameters is discussed more fully below under the heading \emph{\dQuote{Categorical variables}}.
  }
  \item{postrun}{
  Contains a quoted string specifying a command to run (e.g. a linear test) after the model is estimated
  }
  \item{uniform}{
  If uniform=TRUE, random values from a uniform distribution \eqn{unif(-\emph{x}/2,\emph{x}/2)} are added to the perturb variables, where \emph{x} is the \emph{prange} value corresponding with each \emph{pvars} variable. The default is to add values from a normal distribution N(0,\emph{x}.
  }
  \item{niter}{Indicates the number of times to re-estimate the model. Default is 100.}
  \item{object}{a \code{perturb} object to be summarized}
  \item{x}{a \code{summary.perturb} object to be printed}
  \item{dec.places}{number of decimal places to use when printing}
  \item{full}{if TRUE, some extra information is printed}
  \item{\dots}{arguments to be passed on to or from other methods. Print options for class \code{matrix} may be used, e.g. \code{print.gap}}
}

\details{
Perturb is a tool for assessing the impact of small random changes (perturbations) to variables on parameter estimates. It is an alternative for collinearity diagnostics such as \code{\link[car]{vif}} in the \pkg{car} package, \code{\link[rms]{vif}} in the \pkg{rms} package or \code{\link{colldiag}} in this package. Perturb is particularly useful for evaluating collinearity if interactions are present or nonlinear transformations of variables, e.g. a squared term. Perturb can show how the perturbations affect the estimates of a variable and terms derived from it whereas other collinearity diagnostics simply treat interactions and transformations as regular independent variables. Perturb is not limited to linear regression but can be used for all regression-like models. Perturb can also deal with categorical variables by randomly misclassifying them to assess the impact on parameter estimates.

Perturb works by adding a small random \dQuote{perturbation} value to selected independent variables, then re-estimating the model. This process is repeated \emph{niter} times, after which a summary of the means, standard deviation, minimum and maximum of the parameter estimates is displayed. If collinearity is a serious problem in the data, then the estimates will be unstable and vary strongly.

Perturb can be used with categorical variables. Categorical variables are reclassified according to a table of reclassification probabilities. Random numbers determine to which category each case is reclassified at each iteration. The reclassification probabilities are specified to make reclassification to the same category highly likely. For ordered variables, short distance reclassification can be made more likely than long distance. See the section on \emph{\dQuote{Categorical variables}} and \code{\link{reclassify}} for further details.
}

\value{
  An object of class \dQuote{perturb}. The main result is contained in \code{coef.table}, which contains the parameter estimates for each iteration of the perturbation analysis. \code{summary} prints the mean, standard deviation, minimum and maximum of \code{coef.table} over the iterations. If the option \code{full} is specified, \code{\link{reclassify}} prints additional information on how the reclassification probabilities were derived. \code{Summary} also prints information on the transformed model formula.

  \item{coeff.table}{Estimated coefficients for each iteration of the perturbation analysis}
  \item{formula}{The model formula used}
  \item{pvars}{The continuous variables perturbed in the analysis}
  \item{prange}{Magnitude of the perturbations}
  \item{ptrans2}{The transformations using the temporary variables of the perturbation analysis}
  \item{reclassify.tables}{objects produced by \code{\link{reclassify}} for each factor in the perturbation analysis}
  \item{postrun}{A list of data frames with the results of the \code{postrun} specificaton}
  \item{postrun$call}{The \code{postrun} specification}
  \item{formula2}{The model formula using temporary variables}
  \item{distribution}{\dQuote{normal} or \dQuote{uniform}, the distribution from which to draw random numbers}
  \item{perturbed}{A list of data frames with variables that have been modified at each iteration.
  This can be used to recreate the results at a particular iteration}
  \item{message}{\emph{\code{lme4:lmer} only} Warning messages generated when estimating the model}
}

\references{
D. Belsley, E. Kuh, and R. Welsch (1980). \cite{Regression Diagnostics}.
Wiley.

Belsley, D.A. (1991).
\cite{Conditioning diagnostics, collinearity and weak data in regression}. New York: John Wiley & Sons.

Hendrickx, John, Ben Pelzer. (2004). \cite{Collinearity involving ordered and unordered categorical variables}. Paper presented at the RC33 conference in Amsterdam, August 17-20 2004. Available at \url{https://www.researchgate.net/publication/235994590_Collinearity_involving_ordered_and_unordered_categorical_variables}

Hendrickx, John. (2018). \cite{Collinearity in Mixed Models}.
Paper AS03 presented at the \emph{PhUSE EU Connect 2018} conference in Frankfurt, November 4-7 2018.
Available at \url{https://www.lexjansen.com/phuse/2018/as/AS03.pdf}

}

\author{John Hendrickx \email{John.Hendrickx@Danone.com}}

\note{
	Perturb can be used with estimation procedures other than \code{\link{lm}}. On the other hand, collinearity is a result of extreme (multiple) correlation among independent variables. Another option is \code{\link[rms]{vif}} in the \pkg{rms} package \code{\link[car]{vif}} in the \pkg{car} package, or \code{\link{colldiag}}, which use only the independent variables of a regression model. This will usually be a faster solution since maximum likelihood procedures require iterative estimation for each iteration of the perturbation analysis. It is possible though that certain procedures are more sensitive to collinearity than \code{\link{lm}}, in which case perturb could be a better solution.
}

\section{Categorical variables}{
In a perturbation analysis, categorical variables are reclassified with a high probability of remaining in the same category. This could be accomplished in several ways. \code{\link{reclassify}} lets you specify the target percentage reclassification, then adjusts the reclassification frequencies so that the expected frequencies of the reclassified variable are the same as those of the original. In addition, reclassify imposes a meaningful pattern of association between the original and the reclassified variable. See \code{\link{reclassify}} for details.

Categorical variables are specified in \code{perturb} as a list in the \code{pfac} option. The first (unnamed) component is the factor to be reclassified. The names of following components should be valid reclassify options follows by appropriate values. For example, to reclassify the factor \dQuote{type} with a 95\% probability of reclassify to the same category, use:

\code{p2<-perturb(m2,pvars=c("income","education"),prange=c(1,1),}\cr
\code{pfac=list("type",pcnt=95))}

This command will iteratively re-estimate model \code{m2} 100 times (default). Random values from a normal distribution with a standard deviation of 1 will be added to the variables income and education. \code{Reclassify} creates a table of initial reclassification probabilities for \code{type} with a 95\% probability of reclassification to the same category. This initial table is adjusted and the final reclassification probabilities printed in the output are subsequently used to reclassify \code{type} at each iteration.

Use a list of lists to specify a model with more than one reclassification factor, for example:

\code{pfac=list(list("fegp6",pcnt=95),list("eyrc",pcnt=m1),list("expc",pcnt=m2))}\cr
\code{q<-perturb(r1,pfac=pfac)}
}

\section{Post-run commands (e.g. linear tests)}{
The \code{postrun} option can be used to run a command each time the model is re-estimated in \code{perturb}. One reason to do this could be to run linear tests (\emph{\dQuote{contrasts}} in SAS parlance). Take for example the following test from the help section for \code{[car]}\code{\link[car]{linearHypothesis}}:

\code{m1<-lm(prestige~income+education,data=Duncan)}\cr
\code{p3 <- perturb(m1,pvars=c("income","education"),prange=c(1,1),postrun="linearHypothesis(m1, \"income = education\")")}

This will add a \code{postrun} property to \code{p3}. The \code{postrun} can then be extracted, converted to a data frame and values analysed.

\code{# Extract the postrun values}\cr
\code{prun <- p3$postrun}\cr
\code{# Remove the call}\cr
\code{prun$call <- NULL}\cr
\code{# Create a data frame with postrun results}\cr
\code{prun <- do.call("rbind",prun)}

\code{linearHypothesis} creates a data frame with two rows per test, one for the original model, one for the linear test results. The first row is not of interest and can be removed, making use of the fact that all columns except the first two are empty.

\code{# Remove empty rows for the original model}\cr
\code{prun <- prun[!is.na(prun$Df),]}

The results can now be analysed, e.g. by plotting the pvalues:

\code{plot(prun[,6],ylab=colnames(prun)[6])}

Any R function that can be executed on the model object can be used in the postrun option, for example \code{[lmerTest]}\code{\link[lmerTest]{contest}}

\code{q<-perturb(m1,pvars=c("week"),prange=c(0.5),ptrans="week_16<-pmax(week-16,0)",postrun="lmerTest::contest(m1,L)")}

}

\section{Obtaining Model Coefficients}{
The \code{perturb} function uses a \code{coefs} S3 method to obtain model coefficients of models from the packages
\code{multinom}, \code{lme4}, \code{nlme}. These functions may be useful for other purposes as well.

The \code{coefs} function returns a vector of model coefficients. For models estimated by \code{multinom},
the coefficients have been transformed from a matrix to a vector. For models estimated by \code{lme4} or \code{nlme},
the vector contains both fixed and random coefficients.
}

\section{Random Coefficients}{
The \code{VarCorr2} S3 function obtains the random coefficients generated by \code{lme4} and \code{nlme} models
}

\section{Repeated Measures}{
The \code{Repeatedcoefs} S3 function obtains coefficients from repeated measures models by \code{nlme::lme}
}

\seealso{
\code{\link{reclassify}}, \code{\link{colldiag}}, \code{[car]}\code{\link[car]{vif}}, \code{[rms]}\code{\link[rms]{vif}}
}
\examples{
library(car)
data(Duncan)
m1<-lm(prestige~income+education,data=Duncan)
summary(m1)
anova(m1)
vif(m1)
p1<-perturb(m1,pvars=c("income","education"),prange=c(1,1))
summary(p1)
m2<-lm(prestige~income+education+type,data=Duncan)
summary(m2)
anova(m2)
vif(m2)
p2<-perturb(m2,pvars=c("income","education"),prange=c(1,1),pfac=list("type",pcnt=95))
summary(p2)

# Example using lme4:lmer
library(lme4)
data(cd4)
m1<-lmer(logcd4 ~ week+week_16+trt:week+trt:week_16 + (1+week+week_16 | id),data=cd4)
(m1summary<-summary(m1))

# Restrict to 5 iterations for time's sake
q<-perturb(m1,pvars=c("week"),prange=c(3),ptrans="week_16<-pmax(week-16,0)",niter=5)
(q_sum<-summary(q))
# Show warnings generated by models
table(q$messages)

\dontrun{
r1<-lm(ses~fegp6+educyr+eyr+exp2)
summary(r1)
q<-perturb(r1,c("eyr","exp"),c(2.5,2.5),ptrans="exp2<-exp^2")
summary(q)

fegp6<-as.factor(fegp6)

# eyr and exp also as factors
eyrc<-cut(eyr,c(min(eyr),40,50,60,70,80,max(eyr)),include.lowest=T,right=F)
table(eyrc)
expc<-cut(exp,c(0,10,20,max(exp)),include.lowest=T,right=F)
table(expc)

# rough initial reclassification probabilities,
# program will ensure they sum to 100 row-wise
m1<-matrix(0,nlevels(eyrc),nlevels(eyrc))
m1[row(m1)==col(m1)]<-80
m1[abs(row(m1)-col(m1))==1]<-8
m1[abs(row(m1)-col(m1))==2]<-2
m1

m2<-matrix(0,nlevels(expc),nlevels(expc))
m2[row(m2)==col(m2)]<-80
m2[abs(row(m2)-col(m2))==1]<-10
m2[abs(row(m2)-col(m2))==2]<-2
m2

r2<-lm(ses~fegp6+eyrc+expc)
summary(r2)
pfac=list(list("fegp6",pcnt=95),list("eyrc",pcnt=m1),list("expc",pcnt=m2))
q2<-perturb(r2,pfac=pfac,niter=1)
summary(q2)
}
}
\keyword{regression}% at least one, from doc/KEYWORDS
