\name{augSIMEX}
\alias{augSIMEX}
\title{
Analysis of Data with Mixed Measurement Error and Misclassification in Covariates
}
\description{
Implementation of the SIMEX algorithm for data with mixed measurement error and misclassification in covariates.
}
\usage{
augSIMEX(mainformula = formula(data), pimodel = NULL, qimodel = NULL, meformula = NULL, 
        family = gaussian, data, validationdata,
        err.var, mis.var, err.true, mis.true, err.mat = NULL, cppmethod = TRUE,
        repeated = FALSE, repind = list(),
        subset,  offset, weights, na.action, scorefunction=NULL,
        lambda = NULL, M = 50, B = 200, nBoot = 50,
        extrapolation = c("quadratic","linear"),...)
}
\arguments{
  \item{mainformula}{
an object of class ``formula'': an object of class ``formula'': a symbolic description of the model of the response variable, error-prone covariates, and other covariates.
}
  \item{pimodel}{
an object of class ``formula''. A symbolic description in modeling the misclassification rate when binary variable equals to 1. The default choice includes all true covariates in the validation dataset.
}
  \item{qimodel}{
an object of class ``formula''. A symbolic description in modeling the misclassification rate when binary variable equals to 0. The default choice includes all true covariates in the validation dataset.
}
  \item{meformula}{
a list of object of ``formula''. A list of formulas specifying the measurement error model for each error-prone covariate. The default choice will be classic additive models.
}
  \item{family}{
an object of class ``family'' (same as in \code{\link[stats]{glm}} function in stats package).
}
  \item{data}{
a data frame or a matrix of main data. The variable in the main data includes the response, observed covariates matrix which is subject to measurement error, observed binary covariate vector that is prone to misclassification and may also contain the precisely measured covariates matrix. The default choice includes all covariates that mentioned in the \code{mainformula}.
}
  \item{validationdata}{
a data frame or matrix of validation data. The variable in the model includes an observed covariates matrix of that is subject to measurement error and their corresponding precisely measured covariates, an observed binary covariate vector that is prone to misclassification and their corresponding precisely measured covariates, a covariates matrix of precisely measured covariates. 
}
  \item{err.var}{
a vector of character specifying the name of covariates that are subject to measurement error.
}
  \item{mis.var}{
a string specifying the variable name of the binary variable subject to misclassification.
}
  \item{err.true}{
a vector of character specifying the names of the all the counterparts of \code{err.var} that are precisely measure covariates. 
}
  \item{mis.true}{
a string specifying the variable name of the precisely measured binary variable.
}
  \item{err.mat}{
a matrix indicating the variance-covariance matrix of generated Normal distribution in simulation step. 
}
  \item{cppmethod}{
a logical value indicating whether solving the score function via C++ function. The function involves \code{Rcpp} package. The C++ based method is much faster and more computationally efficient than the R based method. See discussion section.
}
  \item{subset}{
a logical vector indicating which subjects should be included in the fitting.
}
  \item{offset}{
a numeric vector indicating the offset. The default choice is null.
}
  \item{weights}{
an optional numeric vector of the weights should be used in the fitting.
}
  \item{na.action}{
a function indicating what method should be applied to deal with missing data.
}
  \item{repeated}{
a logic value indicating whether repeated measurements are involved.
  }
  \item{repind}{
a list of vectors of repeated measurement for error-prone covariates. See detail below.
  }
  \item{scorefunction}{
a function of score function. To allow for the generality, the users can specify their own score function. It should be a function of parameters, response, covariates, weights and offset. A matrix of score value for each individual and each parameter should be returned in the function. An example is shown in \code{\link{glmscore}}
}
  \item{lambda}{
if \code{M} value is not specified, the user can specify a positive sequence of lambda directly. The first element is ususally set to be 0.
}
  \item{M}{
the number of predetermined lambda vector if \code{lambda} is missing. The default is chosen to be 50.
}
  \item{B}{
the number of the dataset generated in the simulation step. The default value is set to be 200.
}
  \item{nBoot}{
the number of the iterations repeated in bootstrap. The default value is set to be 100.
}
  \item{extrapolation}{
specifies the regression model that involves in the extrapolation step. The options include ``linear'' and ``quadratic''. The default is set to be ``quadratic''.
}
  \item{\dots}{
other arguments that pass into the function.
}
}
\details{
In the case of repeated measurements, the users can pick the one of the measurements into the formula of the main model. If more than one covariates have multiple replicates, the users should name the vector of repeated measurements in the list of \code{repind} by the correspoding representitive measurement in the formula of main model. See example.
The number of boostrap should be at least 2. i.e. \code{nBoot}>2.

The examples are mainly for illustration purpose. NA's are possibly to be generated. To obtain precise results, we need to set simulation parameters on a larger scale, which also involves more time for computing. 
}
\value{
\item{coefficient}{the coefficient of the main model after correction.}
\item{se}{a vector of adjusted standard error estimated by bootstrap.}
\item{lambda}{a vector of values that involved in extrapolation step.}
\item{coefmatrix}{a matrix of coefficient before doing extrapolation. Each row correspond to a specified lambda. Each column represent a component of coefficient. This output if for the plotting purpose.}
\item{deviance}{minus twice the maximized log-likelihood.}
\item{null.deviance}{the deviance of null model which only includes intercept and offset term. It is comparable to the deviance. }
\item{df.residual}{the degree of freedom.}
\item{df.null}{the degree of freedom of null model.}
\item{aic}{Akaike's Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters.}
\item{...}{Other components are the arguments that have been used in the function.}


}
\references{
Yi G Y, Ma Y, Spiegelman D, et al. Functional and structural methods with mixed measurement error and misclassification in covariates[J]. Journal of the American Statistical Association, 2015, 110(510): 681-696.
}
\author{
Qihuang Zhang and Grace Y. Yi
}

\seealso{
\code{\link[stats]{glm}}
}
\examples{
### Example 1: Univariate Case
data(ToyUni)
example1<-augSIMEX(mainformula = Y~Xstar+Zstar+W, family = binomial(link=logit),
         data=ToyUni$Main,validationdata=ToyUni$Validation, subset=NULL,
         err.var="Xstar", mis.var="Zstar", err.true="X", mis.true="Z", err.mat = NULL,
         lambda=NULL, M=5, B=2, nBoot=2, extrapolation="quadratic")
summary(example1)

## Without adjustment
example1_naive <- glm(Y~Xstar+Zstar+W, family = binomial(link=logit), data=ToyUni$Main)
summary(example1_naive)

## using  accurate data
example1_true <- glm(Y~Xstar+Zstar+W, family = binomial(link=logit), data=ToyUni$True)


### Example 2: Multivariate Case
data(ToyMult)
meformula<-list(Xstar.X1~-1+X.X1,Xstar.X2~-1+X.X2)    ## measurement error model
example2<-augSIMEX(mainformula = Y~Xstar.X1+Xstar.X2+Zstar+W.W1+W.W2, 
               pimodel=pi~X.X1+X.X2+W.W1+W.W2, qimodel=qi~X.X1+X.X2+W.W1+W.W2, family = binomial,
               data=ToyMult$Main,meformula=meformula,
               validationdata=ToyMult$Validation, subset=NULL,
               err.var=c("Xstar.X1","Xstar.X2"), mis.var="Zstar", err.true=c("X.X1","X.X2"), 
               mis.true="Z", err.mat = NULL,
               lambda=NULL, M=5, B=2, nBoot=2, extrapolation="quadratic")
summary(example2)



### Example 3
data(ToyRepeat)
example3<-augSIMEX(mainformula = Y~Xstar1+Zstar+W, family = binomial(link=logit),
            data=ToyRepeat$Main,validationdata=ToyRepeat$Validation, 
            subset=NULL, err.var="Xstar1", mis.var="Zstar", err.true="X", mis.true="Z", 
            err.mat = NULL, repeated = TRUE,repind=list(Xstar1=c("Xstar1","Xstar2")),
            lambda=NULL, M=5, B=2, nBoot=2, extrapolation="quadratic")
summary(example3)
}
\keyword{funtion}
\keyword{core}
