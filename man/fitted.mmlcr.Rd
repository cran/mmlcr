\name{fitted.mmlcr}
\alias{fitted.mmlcr}
\alias{residuals.mmlcr}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Extract mmlcr Model Component Fitted or Residual Values }
\description{
  Extract fitted or residual values from components of a mmlcr object.
}
\usage{
fitted.mmlcr(object, which = 1:length(object$components), ...)
residuals.mmlcr(object, which = 1:length(object$components), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ Fitted mmlcr object. }
  \item{which}{ Which components to extract fit or residuals from}
  \item{\dots}{ Included for compatibility with generic functions. }
}
\details{
  For each observation in a component, the fit is the weighted mean of
  the fits from the individual class fits. The weights are the posterior class
  membership probabilities.
}
\value{
  A list, with a vector of fits or residuals for each component in which.
}

\section{WARNING }{Although the function works on multinom components,
the results are not meaningful.}

\keyword{models}

