\name{summary.augSIMEX}
\alias{summary.augSIMEX}

\title{
Summarizing the Adjusted Fits of Generalized Linear Model 
}
\description{
This function is a method for \code{augSIMEX} object.
}
\usage{
## S3 method for class 'augSIMEX'

\method{summary}{augSIMEX}(object, dispersion = NULL, ...)
}

\arguments{
  \item{object}{
  the ``augSIMEX'' object gotten from \code{\link[augSIMEX]{augSIMEX}} function.
}
  \item{dispersion}{
  the dispersion parameter for the family used.
}
  \item{\dots}{
  other arguments passed to the function.
}
}
\details{
The function provides summary output similar to \code{\link[stats]{glm}}. 
}
\value{
\item{coefficients}{the matrix of coefficients, standard errors, z-values and p-values.}
\item{dispersion}{the component from arguments of the function}
\item{call}{the component from \code{object}}
\item{family}{the component from \code{object}}
\item{deviance}{the component from \code{object}}
\item{aic}{the component from \code{object}}
\item{df.residual}{the component from \code{object}}
\item{null.deviance}{the component from \code{object}}
\item{df.null}{the component from \code{object}}
\item{extrapolation}{the component from \code{object}}
\item{err.var}{the component from \code{object}}
\item{mis.var}{the component from \code{object}}
\item{err.true}{the component from \code{object}}
\item{mis.true}{the component from \code{object}}
\item{B}{the component from \code{object}}
\item{mis.true}{the component from \code{object}}
\item{nBoot}{the component from \code{object}}
}

\author{
 Qihuang Zhang and Grace Y. Yi.
}

\seealso{
\code{\link[stats]{summary.glm}},\code{\link[augSIMEX]{augSIMEX}}
}

\examples{
## Please see the example in augSIMEX function
}
\keyword{function}

