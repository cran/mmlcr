\name{plot.mmlcronce}
\alias{plot.mmlcronce}

\title{Plot a mmlcronce Component}
\description{
  Generates a plot for mmlcronce component of a mmlcr object, showing weighted means plus and minus one standard deviation for
  each class.
   Generally called by \code{plot.mmlcr}.
}
\usage{
plot.mmlcronce(x, post.prob, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ A mmlcr object. }
  \item{post.prob}{ Data.frame of class membership posterior probabilities. }
  \item{\dots}{ Other arguments to be passed to \code{plot}.  }
}
\note{
 There is no component plot for the multinomonce components.
}

\seealso{ \code{\link{plot.mmlcr}}, \code{\link{mmlcr}}, \code{\link{mmlcrObject}} }
\keyword{models}
\keyword{hplot}
