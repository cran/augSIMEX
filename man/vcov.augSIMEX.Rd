\name{vcov.augSIMEX}
\alias{vcov.augSIMEX}

\title{
Extract the variance-covariance matrix from the fitted augSIMEX object
}
\description{
This function return the variance-covariance matrix from a fitted augSIMEX object
}
\usage{
## S3 method for class 'augSIMEX'

\method{vcov}{augSIMEX}(object, ...)
}

\arguments{
  \item{object}{
the ``augSIMEX" object gotten from \code{\link[augSIMEX]{augSIMEX}} function.
}
  \item{...}{ other arguments to pass to the function.
}
}

\author{
Qihuang Zhang and Grace Y. Yi.
}

\seealso{
\code{\link[stats]{vcov}}
}

\examples{
data(ToyUni)
example <- augSIMEX(mainformula = Y ~ Xstar + Zstar + W, family = binomial(link = logit),
  mismodel = pi|qi ~ W, 
  meformula = Xstar ~ X + Z + W,
  data = ToyUni$Main,validationdata = ToyUni$Validation, subset = NULL,
  err.var = "Xstar", mis.var = "Zstar", err.true = "X", mis.true = "Z", 
  err.mat = NULL,
  lambda = NULL, M = 5, B = 2, nBoot = 2, extrapolation="quadratic")
vcov(example)
}

\keyword{function}
