
\name{BIC.mmlcr}
\title{Bayesian Information Criterion}
\usage{
\method{BIC}{mmlcr}(object, ...)
}
\alias{BIC.mmlcr}

\arguments{
 \item{object}{a fitted mmlcr object.}
 \item{\dots}{optional fitted model objects.}
}
\description{
  This generic function calculates the Bayesian information criterion,
  also known as Schwarz's Bayesian criterion (SBC), for an mmlcr object for which a log-likelihood value can be obtained,
  according to the formula \eqn{-2 \mbox{log-likelihood} + n_{par}
    \log(n_{obs})}{-2*log-likelihood + npar*log(nobs)}, where
  \eqn{n_{par}}{npar}  represents the
  number of parameters and \eqn{n_{obs}}{nobs} the number of
  observations in the fitted model. 
}
\value{
 if just one object is provided, returns a numeric value with the
 corresponding BIC; if more than one object are provided, returns a
 \code{data.frame} with rows corresponding to the objects and columns
 representing the number of parameters in the model (\code{df}) and the
 BIC.
}
\references{
  Schwarz, G. (1978) "Estimating the Dimension of a Model", Annals of
  Statistics, 6, 461-464.
}
\seealso{\code{\link{AIC}}, \code{\link{mmlcrObject}}}
\examples{
\dontrun{data(mmlcrdf)}
\dontrun{mmlcrdf.mmlcr2 <- mmlcr(outer =  ~ sex + cov1 | id, 
components = list(
	list(formula = resp1 ~ 1, class = "cnormonce", min = 0, max = 50),
	list(formula = resp2 ~ poly(age, 2) + tcov1, class = "poislong"),
	list(formula = resp3 ~ poly(age, 2), class = "multinomlong")
), data = mmlcrdf, n.groups = 2)}

\dontrun{BIC(mmlcrdf.mmlcr2)}
}
\keyword{models}

