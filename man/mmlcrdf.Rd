\name{mmlcrdf}
\alias{mmlcrdf}
\title{
An artificial data frame for `mmlcr' examples.
}
\description{
The \code{mmlcrdf} data frame has 750 rows and 8 columns. It is a simulated data frame
used for the \code{mmlcr} examples.
}
\usage{
data(mmlcrdf)
}
\format{
\describe{
This data frame contains the following columns:
\item{\code{id}}{subject id number. Each subject has 5 observations.}
\item{\code{age}}{subject's age}
\item{\code{sex}}{subject's sex}
\item{\code{cov1}}{a continuous covariate, fixed over time}
\item{\code{tcov1}}{a time-varying, continuous covariate}
\item{\code{resp1}}{a continuous response variable, left-censored at 0}
\item{\code{resp2}}{a count variable}
\item{\code{resp3}}{a factor response variable}
}}
\source{
simulation
}
\seealso{
\code{\link{mmlcr}}
}
\examples{
\dontrun{data(mmlcrdf)}
\dontrun{mmlcrdf.mmlcr2 <- mmlcr(outer =  ~ sex + cov1 | id, 
	components = list(
		list(formula = resp1 ~ 1, class = "cnormonce", min = 0, max = 50),
		list(formula = resp2 ~ poly(age, 2) + tcov1, class = "poislong"),
		list(formula = resp3 ~ poly(age, 2), class = "multinomlong")
	), data = mmlcrdf, n.groups = 2)}
}
\keyword{datasets}
