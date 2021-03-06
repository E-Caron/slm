% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slm-method.R
\name{summary.slm}
\alias{summary.slm}
\title{Summarizing Stationary Linear Model Fits}
\usage{
\method{summary}{slm}(object, correlation = FALSE,
  symbolic.cor = FALSE, ...)
}
\arguments{
\item{object}{an object of class "\code{slm}", usually, a result of a call to \code{slm}.}

\item{correlation}{logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed.}

\item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in a symbolic form (see \code{\link[stats:symnum]{symnum}}) rather than as numbers.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The function \code{summary.slm} computes and returns a list of summary statistics of the fitted linear model
 given in \code{object}, using the components (list elements) "\code{call}" and "\code{terms}" from its argument, plus:
 \item{residuals }{ the residuals, that is response minus fitted values.}
 \item{coefficients }{ a \eqn{p*4} matrix with columns for the estimated coefficient, its standard error, z-statistic and corresponding (two-sided) p-value.
  Aliased coefficients are omitted.}
 \item{aliased }{ named logical vector showing if the original coefficients are aliased.}
 \item{sigma }{ the square root of the estimated variance of the error process.}
 \item{df }{ degrees of freedom, a 3-vector \eqn{(p, n-p, p*)}, the first being the number of non-aliased coefficients, the last being the total number of coefficients.}
 \item{chi2statistic }{ a 2-vector with the value of the chi2-statistic with its degree of freedom.}
 \item{r.squared }{ \eqn{R^2}, the 'fraction of variance explained by the model'.}
 \item{cov.unscaled }{ the matrix (X^{t} X)^{(-1)}.}
 \item{correlation }{ the correlation matrix corresponding to the above \code{cov.unscaled}, if \code{correlation = TRUE} is specified.}
 \item{symbolic.cor }{ (only if \code{correlation} is true.) The value of the argument \code{symbolic.cor}.}
}
\description{
Summary method for class "\code{slm}".
}
\examples{
data("shan")
reg1 = slm(shan$PM_Xuhui ~ . , data = shan, method_cov_st = "fitAR", model_selec = -1)
summary(reg1)

data("co2")
y = as.vector(co2)
x = as.vector(time(co2)) - 1958
reg2 = slm(y ~ x + I(x^2) + I(x^3) + sin(2*pi*x) + cos(2*pi*x) + sin(4*pi*x) +
 cos(4*pi*x) + sin(6*pi*x) + cos(6*pi*x) + sin(8*pi*x) + cos(8*pi*x),
 method_cov_st = "fitAR", model_selec = -1)
summary(reg2)
}
\references{
E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
 \url{https://arxiv.org/abs/1906.06583}.
}
\seealso{
The model fitting function \code{\link[slm:slm]{slm}}, \code{\link[base:summary]{summary}}.

 The function \code{\link[stats:coef]{coef}} extracts the matrix of coefficients with standard errors, z-statistics and p-values.
}
