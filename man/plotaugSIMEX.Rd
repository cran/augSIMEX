\name{plotaugSIMEX}
\alias{plotaugSIMEX}

\title{
Plot of Extrapolation
}
\description{
This function visualizes the extrapolation step. It plots the simulated value for each lambda value and the curve of extrapolation. The function can provide the plots for both methods simultaneously, and guides the user to choose proper extrapolation method. 
}
\usage{
## S3 method for class 'augSIMEX'

plotaugSIMEX(object, variable, extrapolation = c("both", "linear", "quadratic"),
       xlim = c(-1.2, 1.2), ylim = c(-3, 3), ...)
}

\arguments{
  \item{object}{
the ``augSIMEX" object gotten from \code{\link[augSIMEX]{augSIMEX}} function.
}
  \item{variable}{
 a value or a vector specifying which variable the plot is referred to. It can be the index of the variable in the dataset to be plotted, or it can be the name of the true variable (specified in err.ture).
}
  \item{extrapolation}{
 the extrapolation methods used in the plot choosing from  ``quadratic'', ``linear'' or ``both''.
}
  \item{xlim}{
the range of x-axis in the plot.
}
  \item{ylim}{
the range of y-axis in the plot.
}
  \item{\dots}{
other arguments that are passed into the function.
}
}
\details{
The user may need to adjust the range of y axis for the proper display.
}
\author{
Qihuang Zhang and Grace Y. Yi.
}

\seealso{
\code{\link{plot}}
}

\examples{
data(ToyUni)
example<-augSIMEX(mainformula = Y~Xstar+Zstar+W, family = binomial(link=probit),
         data=ToyUni$Main,validationdata=ToyUni$Validation, subset=NULL,
         err.var="Xstar", mis.var="Zstar", err.true="X", mis.true="Z", err.mat = NULL,
         lambda=NULL, M=5, B=2, nBoot=2, extrapolation="quadratic")
         
plotaugSIMEX(example,"X",extrapolation = "both",ylim = c(-1,0.4) )
}

\keyword{function}
