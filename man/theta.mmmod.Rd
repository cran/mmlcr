\name{theta.mmmod}
\alias{theta.mmmod}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Estimate theta of the Negative Binomial by Moments }
\description{
  Estimate theta of the negative binomial by moments, adapted to
  mmlcrfit from theta.mm of the MASS library
}
\usage{
theta.mmmod(y, u, dfr, wts, limit = 10, eps = .Machine$double.eps^0.25)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ Vector of observed values. }
  \item{u}{ Estimated mean vector. }
  \item{dfr}{ Residual degrees of freedom. }
  \item{wts}{ Weights for the observations. }
  \item{limit}{ Limit on iterations. }
  \item{eps}{ Convergence tolerance. }
}

}
\value{
  Estimate of theta.
}

\seealso{ \code{\link{theta.mm}},  \code{\link{mmlcrfit}}, \code{theta.mm} from the MASS library}

\keyword{distribution}
\keyword{models}

