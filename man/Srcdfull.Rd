\name{Srcdfull}
\alias{Srcdfull}
\title{ Longitudinal Dataset of Aggression and Reading Ability in Children. }
\usage{data(Srcdfull)}
\description{
The \code{Srcdfull} data frame has 1620 rows and 9 columns.
It is a longitudinal dataset from 405 children  with scores for aggression 
and reading ability.}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{Subject id.}
    \item{sex}{A factor variable with the sex of each subject.}
    \item{age}{Age of subject at time of interview.}
    \item{anti}{Antisocial score, an integer from 0 to 12.}
    \item{read}{Reading ability score, ranging from 0 to 8.4.}
    \item{momage}{Mother's age at first interview.}
    \item{kidage}{Subject's age at first interview.}
    \item{homecog}{An integer score, ranging from 0 to 14, 
        at the first interview for cognitive stimulation at home.}
    \item{homeemo}{An integer score, ranging from 0 to 14, 
        at the first interview for emotional support at home.}
  }
}
\details{
 The sample consists of N=405 children drawn from the Children of 
 the National Longitudinal Survey of Youth, about half of whom are 
 missing one or more of the repeated measures on aggression or reading 
 ability. 
 The measures of interest include four repeated measures of aggressive 
 behavior taken at two-year intervals, four repeat measures of reading 
 achievement taken at two-year intervals, initial interview measures of 
 cognitive stimulation and social support of the child in the home, and 
 background measures of child gender, child age, and mother age.
}
\source{
\url{http://www.unc.edu/~curran/srcd-docs/srcddata.txt}
}
\references{
\url{http://www.unc.edu/~curran/srcd-docs/srcdmeth.pdf}
}
\examples{
data(Srcdfull)
Srcdfull.mmlcr4 <- mmlcr(outer = ~1 | id, components = list(
    list(formula = anti ~ poly(age, 2), min = 0, max = 12, 
    class = "cnormlong"), 
    list(formula = read ~ poly(age, 2), class = "normlong")), 
    data = Srcdfull, n.groups = 4)
}
\keyword{datasets}

