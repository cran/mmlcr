\name{plot.mmlcr}
\alias{plot.mmlcr}
\alias{plot.multinomlong}
\alias{plot.multinomonce}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plot Components of mmlcr Object}
\description{
  Generates a plot for each component of a mmlcr object.
  Longitudinal components are plotted against time for each class, 
  while ``once'' components show means and standard deviations for
  each class.
}
\usage{
plot.mmlcr(x, which = 1:length(x$components), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ A mmlcr object. }
  \item{which}{ Which components to plot (by number). }
  \item{\dots}{ Arguments to be passed to component plots. In particular, a ``smooth'' parameter, equal to \code{supsmu}'s \code{bass} can be used in 
longitudinal components.  }
}
\details{
  Longitudinal components show plots of means for each class over the
  first covariate.  Non-longitudinal components show the weighted mean plus and minus one weighted standard deviation for each class. There is no component plot for the multinom components.
}

\seealso{ \code{\link{mmlcr}}, \code{\link{plot.mmlcrlong}}, \code{\link{plot.mmlcronce}}, \code{\link{mmlcrObject}} }
\keyword{models}
\keyword{hplot}
