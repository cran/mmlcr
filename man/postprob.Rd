\name{postprob}
\alias{postprob}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Extract Posterior Probabilities from mmlcr Object}
\description{
  Extract posterior class probabilities from an mmlcr object.
}
\usage{
postprob(object)
}
\arguments{
  \item{object}{ A mmlcr object. }
}
\value{
  A data.frame of posterior probabilities of class membership
}

\seealso{  \code{\link{mmlcr}}, \code{\link{mmlcrObject}} }


\examples{
\dontrun{
data(mmlcrdf)

mmlcrdf.mmlcr2 <- mmlcr(outer =  ~ sex + cov1 | id, 
components = list(
	list(formula = resp1 ~ 1, class = "cnormonce", min = 0, max = 50),
	list(formula = resp2 ~ poly(age, 2) + tcov1, class = "poislong"),
	list(formula = resp3 ~ poly(age, 2), class = "multinomlong")
), data = mmlcrdf, n.groups = 2)

postprob(mmlcrdf.mmlcr2)

}}
\keyword{models}
