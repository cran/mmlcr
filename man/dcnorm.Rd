\name{dcnorm}
\alias{dcnorm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Density function for Censored Normal Distribution }
\description{
  Density function for the censored (not truncated!) normal distribution
  in a form suitable for mmlcrlike
}
\usage{
dcnorm(x, mean = rep(0, length(x)), sigma = rep(1,
                 length(x)), min = rep(-(.Machine$double.xmax),
                 length(x)), max = rep(.Machine$double.xmax,
                 length(x)))
}
\arguments{
  \item{x}{ vector of quantiles }
  \item{mean}{ vector of means}
  \item{sigma}{ vector of standard deviations }
  \item{min}{ vector of left censoring cutoffs }
  \item{max}{ vector of right censoring cutoffs }
}

\value{
  A vector of densities.
}

\seealso{\code{\link{mmlcr}}, \code{\link{mmlcrlike}} }

\examples{

dcnorm(0, mean = 0, sigma = 1)

dcnorm(0, mean = 0, sigma = 1, min = 0)
}
\keyword{distribution}
