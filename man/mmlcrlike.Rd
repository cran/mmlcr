\name{mmlcrlike}
\alias{mmlcrlike}
\alias{mmlcrlike.cnormlong}
\alias{mmlcrlike.default}
\alias{mmlcrlike.multinomlong}
\alias{mmlcrlike.nb1long}
\alias{mmlcrlike.nblong}
\alias{mmlcrlike.normlong}
\alias{mmlcrlike.poislong}
\title{
Used internally by mmlcr to calculate component likelihoods. 
}
\description{Find the likelihoods for each subject in an \code{mmlcr} component object. Meant to be
internally by the function \code{mmlcr}
}
\usage{
mmlcrlike(object, weights, classnumber, final = FALSE)

}
\arguments{
\item{object}{
a component of an \code{mmlcr} object
}
\item{weights}{
weights for the observations in object
}
\item{classnumber}{
the class number for the likelihoods
}
\item{final}{
see below
}}
\value{

if final = FALSE, a vector of likelihoods corresponding to the
weights. The weights are not used except for their row.names.


if final = TRUE, the component is returned, with fitted and residual
components added. The fitted and residual components will only be
correct if mmlcrlike is called sequentially over the full set of 
classnumbers, as is the case when mmlcrlike is called by \code{mmlcr}.
}
\section{Warning}{
the mmlcrlike.multinomlong method returns fitted and residual components
that have no real meaning. It is included only for completeness.
}
\seealso{
\code{\link{mmlcr}}, \code{\link{mmlcrObject}}
}
\keyword{models}
