\name{dnb2}
\alias{dnb2}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Density for negative binomial, used in mmlcr}
\description{
The density for the negative binomial function, in a form suitable for the mmlcrlike function.}
\usage{
dnb2(y, mu, th)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{Observations }
  \item{mu}{The mean of the distribution }
  \item{th}{ The parameter \eqn{\theta}. }
}
\details{
  The variance function of the distribution is \eqn{\mu +\mu^2/\theta}.
  
  The density is 
  \deqn{
  f(y) = \frac{\Gamma(\theta + y) \mu^y \theta^\theta}{\Gamma(\theta) y! (\mu + \theta)^(\mu+\theta)}.
}{%
f(y) = [Gamma(th + y) mu^y theta^th] 
		/ [Gamma(th) y! (mu + th)^(mu + th)].
}
}
\value{
  The density.
}
\references{Venables, W. N. and Ripley, B. D. 2002 \emph{Modern Applied Statistics with S}, Springer Verlag. }

\seealso{ \code{\link{mmlcr}}, \code{\link{mmlcrlike}} }
\keyword{distribution}
