\name{plot.mmlcrlong}
\alias{plot.mmlcrlong}

\title{Plot a mmlcrlong Component}
\description{
  Generates a plot for mmlcrlong component of a mmlcr object.
  Class means are plotted against ``time.'' Generally called by \code{plot.mmlcr}.
}
\usage{
plot.mmlcrlong(x, post.prob, smooth = 0, xlab = names(data)[2], ylab = 
	names(data)[1], ylim = c(min(yy, na.rm = TRUE), 
	max(yy, na.rm = TRUE)), cols = rep(1, dim(post.prob)[2]), 
	pch = as.character(1:(dim(post.prob)[2])), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ A mmlcr object. }
  \item{post.prob}{ Data.frame of class membership posterior probabilities. }
  \item{smooth}{  A smoothing parameter, equal to \code{supsmu}'s \code{bass}.}
  \item{xlab}{ Label for the x-axis. }
  \item{ylab}{ Label for the y-axis. }
  \item{ylim}{ Limits for the y-axis. }
  \item{cols}{ The cols argument for \code{par}. }
  \item{pch}{ The pch argument for \code{par}. }
  \item{\dots}{ Other arguments to be passed to \code{plot}.  }
}
\note{
 There is no component plot for the multinomlong components.
}

\seealso{ \code{\link{plot.mmlcr}}, \code{\link{mmlcr}}, \code{\link{mmlcrObject}} }
\keyword{models}
\keyword{hplot}
