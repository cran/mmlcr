\name{mmlcrcomponentinit}
\alias{mmlcrcomponentinit}
\alias{mmlcrcomponentinit.cnormlong}
\alias{mmlcrcomponentinit.cnormonce}
\alias{mmlcrcomponentinit.default}
\alias{mmlcrcomponentinit.multinomlong}
\alias{mmlcrcomponentinit.multinomonce}
\alias{mmlcrcomponentinit.nblong}
\alias{mmlcrcomponentinit.nbonce}
\alias{mmlcrcomponentinit.normlong}
\alias{mmlcrcomponentinit.normonce}
\alias{mmlcrcomponentinit.poislong}
\alias{mmlcrcomponentinit.poisonce}    
\alias{summary.cnormlong}
\alias{summary.multinomlong}
\alias{summary.nblong}
\alias{summary.normlong}
\alias{summary.poislong}

\title{
Initialize an `mmlcr' component
}
\description{


Used internally by \code{mmlcr} to initialize an \code{mmlcr} component. The function is generic.
}
\usage{
mmlcrcomponentinit(object, n.groups,  prob, data, grouping)
}
\arguments{
\item{object}{
	see details below
}
\item{n.groups}{
	the number of classes. Provided by \code{mmlcr}.
}
\item{prob}{
	a data.frame of (initial) posterior probabilities of class
	membership. Provided by \code{mmlcr}.
}
\item{data}{
	a data.frame. Provided by \code{mmlcr}.
}
\item{grouping}{
	the variable of data indicating subject id. Provided by \code{mmlcr}.
}}
\value{


an \code{mmlcr} component object.
}
\details{


The object argument is a list, originally provided to the \code{mmlcr} function as a 
component of a list in the components argument. At a minimum, that list consists 
of a formula of the form \code{resp ~ covars}, where resp is the response variable for 
that component and covars are any additional covariates (\code{resp ~ 1} is commonly 
used), and a class. For longitudinal responses, the covariates should include 
the chronological variable, e.g., \code{resp ~ poly(age, 2)}. Implemented classes include 
\code{cnormlong}, \code{cnormonce}, \code{multinomlong}, \code{multinomonce}, \code{nblong}, \code{nbonce}, 
\code{normlong}, \code{normonce}, \code{poislong}, and \code{poisonce}. The *\code{long} classes are for longitudinal
responses, while the *\code{once} classes are for one-time responses. So that a rectangular
data structure can be used for the data argument in \code{mmlcr}, the *\code{once} classes use
just the first occurence of a response for each subject.


The \code{cnorm}* classes are for censored normal distributions. They take a \code{min} and \code{max}
argument for the minimum and maximum possible observed value. The \code{multinom}* classes
are for multinomial, or categorical responses. The response variable should be
a factor.  The \code{nb}* classes are for 
negative binomial models (the usual overdispersed Poisson model). The \code{norm}* classes 
are for normal models, and the \code{pois}* classes are for Poisson models. For the \code{cnorm}* 
and \code{norm}* classes the minimum and maximum variances are constrained to be within a 
factor of 100 of each other. For the \code{pois}* classes an overdispersion estimate is 
calculated but not used.
}

\seealso{




\code{\link{mmlcr}}, \code{\link{mmlcrObject}}
}
\keyword{models}

