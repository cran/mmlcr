\name{mmlcrfit}
\alias{mmlcrfit}
\alias{mmlcrfit.cnormlong}
\alias{mmlcrfit.default}
\alias{mmlcrfit.multinomlong}
\alias{mmlcrfit.nblong}
\alias{mmlcrfit.nofit}
\alias{mmlcrfit.normlong}
\alias{mmlcrfit.poislong}
\title{
Used internally by mmlcr to fit mmlcr component objects.
}
\description{
fits an \code{mmlcr} component object. Used internally by the 
function \code{mmlcr}
}
\usage{
mmlcrfit(object, weights, classnumber)
}
\arguments{
\item{object}{
an \code{mmlcr} component
}
\item{weights}{
a data.frame of weights, with row.names indicating subject ids.
}
\item{classnumber}{
a number indicating the latent class to use, i.e., which column of weights to use.
}}
\value{
the same \code{mmlcr} component object with an updated fit.
}
\seealso{
\code{\link{mmlcr}}, \code{\link{mmlcrObject}}, \code{\link{mmlcrlike}}
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
