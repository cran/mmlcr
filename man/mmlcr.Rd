\name{mmlcr}
\alias{mmlcr}
\alias{mmlcr.default}
\alias{mmlcr.mmlcr}
\alias{anova.mmlcr}
\alias{formula.mmlcr}
\alias{logLik.mmlcr}
\alias{print.mmlcr}
\alias{print.summary.mmlcr}
\alias{summary.mmlcr}
\alias{vcov.mmlcr}         


\title{
Mixed Mode Latent Class Regression
}
\description{
This function fits a latent class regression model, where the responses may be of
different types, as well as longitudinal or time-invariant.
}
\usage{
mmlcr(object = NULL, ...)
\method{mmlcr}{default}(object = NULL, outer, components, data, 
	subset = TRUE, n.groups = 2, prior.prob = NULL, 
	post.prob = NULL, no.fit = FALSE, max.iter = 50, trace = TRUE, 
	tol = 0.005, ...)
\method{mmlcr}{mmlcr}(object, max.iter = 50, trace = TRUE, tol = 0.005, ...)
}
\arguments{
\item{outer}{
a formula of the form  \code{~ predictors | id} where the model statement for predictors is
as for regression models and \code{id} is the subject factor.  The expression \code{1} can be used 
for the predictors when there is no other suitable candidate. 
}
\item{components}{
a list of the response components. Each element of the list is itself a list. 
At a minimum, that list consists of a formula of the form \code{resp ~ covars}, where \code{resp} is the
response variable for that component and \code{covars} are any additional 
covariates (\code{resp ~ 1} is commonly used, or something like \code{resp ~ age} for longitudinal classes), 
and a class. Implemented classes include 
\code{cnormlong}, \code{cnormonce}, \code{multinomlong}, \code{multinomonce}, \code{nb1long},  \code{nblong}, \code{nbonce}, 
\code{normlong}, \code{normonce}, \code{poislong}, and \code{poisonce}. Some classes require additional list 
elements. Details can be found with the \code{mmlcrcomponentinit} help files.
}
\item{data}{
a data frame in which to interpret the variables occurring in outer and in the 
formulas for the individual components.
}

\item{subset }{
expression giving which subset of the rows of the data should be used in the fit. 
All observations are included by default.
}
\item{n.groups}{
the number of latent classes.
}
\item{prior.prob }{
a data.frame giving the prior probabilities of class membership according 
to covariate information, not ``prior knowledge.'' The row.names should
match the subject id given in outer. If no value is given, which is recommended, the
default is equal probabilities.
}
\item{post.prob }{
a data.frame giving the posterior probabities of class membership. The row.names 
should match the subject id given in outer. If \code{post.prob} is not given, starting values
are assigned randomly. If one is fitting a modification of a prior fit called \code{fit1}, 
perhaps by changing the covariates, the recommendation is to include the term
\code{post.prob = fit1$post.prob}.
}
\item{no.fit }{
if TRUE, returns a \code{mmlcrObject} without fitting. The \code{mmlcrObject}
could be fit later by a call to \code{mmlcr(mmlcrObject)}.
}
\item{object}{
an \code{mmlcr} object. See the \code{mmlcrObject} help file for details.
}
\item{max.iter}{
the maximum number of iterations
}
\item{trace}{
if TRUE, traces through the iterations, giving the loglikelihood, a convergence speed
index, the loglikelihood goal, and the current class percentages.
}
\item{tol}{
the tolerance between the loglikelihood and the loglikelihood goal, used as a
stopping criterion.
}
\item{...}{possible additional arguments.}
}
\value{
a fitted \code{mmlcrObject}
}
\details{
The components portion of the call refers to components of the response, not
components in the sense, used in much of the mixtures literature, that is here
called latent classes or groups.


It is not yet possible to specify a prior, or even starting values, for the parameters
of the individual latent classes. Instead, one can assign starting values to the posterior
probabilities of each individual via the post.prob part of the function call.
}
\seealso{
\code{\link{mmlcrObject}}, \code{\link{mmlcrcomponentinit}}, 
\code{\link{summary.mmlcr}}, \code{\link{plot.mmlcr}}
}
\examples{
data(mmlcrdf)

mmlcrdf.mmlcr2 <- mmlcr(outer =  ~ sex + cov1 | id, 
components = list(
	list(formula = resp1 ~ 1, class = "cnormonce", min = 0, max = 50),
	list(formula = resp2 ~ poly(age, 2) + tcov1, class = "poislong"),
	list(formula = resp3 ~ poly(age, 2), class = "multinomlong")
), data = mmlcrdf, n.groups = 2)


mmlcrdf.mmlcr2.inter <- mmlcr(outer =  ~ sex * cov1 | id, 
components = list(
	list(formula = resp1 ~ 1, class = "cnormonce", min = 0, max = 50),
	list(formula = resp2 ~ poly(age, 2) + tcov1, class = "poislong"),
	list(formula = resp3 ~ poly(age, 2), class = "multinomlong")
), data = mmlcrdf, n.groups = 2,
post.prob = mmlcrdf.mmlcr2$post.prob, no.fit = TRUE)


mmlcrdf.mmlcr2.inter <- mmlcr(mmlcrdf.mmlcr2.inter)
}
\keyword{models}
% Converted by Sd2Rd version 1.21.






