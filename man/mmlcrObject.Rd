\name{mmlcrObject}
\alias{mmlcrObject}
\title{
Mixed Mode Latent Class Regression Object
}
\description{
Classes of objects returned by fitting mixed mode latent class regression objects.  This class of objects is returned by the \code{mmlcr} function to represent a fitted 
mixed mode latent class regression model.  Objects of this class have methods for the functions \code{print}, \code{plot}, \code{summary}, \code{anova}, \code{vcov}, \code{fitted}, and \code{residuals}, among others.
}
\details{
The following components must be included in a legitimate \code{mmlcr} object. 
\item{components}{
a list of the response components. Each element of the list is itself a 
list, which includes the data for that response, the fit coefficients, 
the fitted response and the response residuals, among others.
}
\item{outer}{
a formula for latent class membership.
}
\item{outer.df}{
a data.frame with the variables used in outer.
}
\item{gamma.matrix}{
the coefficients from the multinomial fit for latent classes. The values 
are log(odds ratios) compared to the first class.
}
\item{prior.prob}{
a data.frame giving the prior probabilities (that is, based on the fit from
outer) of class membership for each subject.
}
\item{post.prob}{
a data.frame giving the posterior probabilities (that is, based on the 
prior probabilities as well as the data from the components) of class 
membership for each subject.
}
\item{loglikelihood}{
the loglikelihood of the final fit.
}
\item{BIC}{
the BIC of the final fit.
}
\item{df}{
the degrees of freedom of the model
}
\item{resid.df}{
the residual degrees of freedom of the model
}
\item{call}{
the call to the \code{mmlcr} function
}}
\seealso{
\code{\link{mmlcr}}
}
\keyword{models}
% Converted by Sd2Rd version 1.21.

