\name{mmlcrclassify}
\alias{mmlcrclassify}
\title{
Classification of Individuals by Modal Posterior Probabilities
}
\description{
This function gives a latent class classification for each individual by
giving the class to which the individual has the highest posterior probabilities.
}
\usage{
mmlcrclassify(object)
}
\arguments{
\item{object}{an \code{mmlcrObject}.}
}
\value{
a data.frame with row.names equal to the ids and with a single column giving the
classification
}
\seealso{
\code{\link{mmlcrObject}}, \code{\link{mmlcr}}
}
\examples{
\dontrun{
data(mmlcrdf)

mmlcrdf.mmlcr2 <- mmlcr(outer =  ~ sex + cov1 | id, 
components = list(
	list(formula = resp1 ~ 1, class = "cnormonce", min = 0, max = 50),
	list(formula = resp2 ~ poly(age, 2) + tcov1, class = "poislong"),
	list(formula = resp3 ~ poly(age, 2), class = "multinomlong")
), data = mmlcrdf, n.groups = 2)

mmlcrclassify(mmlcrdf.mmlcr2)}
}
\keyword{models}

