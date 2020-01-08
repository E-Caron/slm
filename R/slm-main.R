#document() for Roxygen documentation
#S3 $ -> attributes
#S4 @ : reg@cov_st
#?slm for documentation
#devtools::build_manual
#or setwd("/Users/emmanuel/Universite-Recherche/Seafile/Work-PackageR-STemp/slm/")
#system("R CMD Rd2pdf . --title=Package slm --output=./manual.pdf --force --no-clean --internals")

#' @title slm class
#'
#' @description An S4 class to create an \code{slm} object.
#'
#' @slot method_cov_st the method used to compute the autocovariance vector of the error process.
#' @slot cov_st a numeric vector with the estimated autocovariances of the error process, computed from
#' the \code{method_cov_st} method.
#' @slot Cov_ST the estimated covariance matrix of the error process, computed from the \code{method_cov_st} method.
#' @slot model_selec the order of the chosen method. If \code{model_selec = -1}, the method works automatically.
#' @slot norm_matrix the normalization matrix of the design X.
#' @slot design_qr the matrix \eqn{(X^{t} X)^{-1}}.
#' @inheritParams lm
#@slot lm contains the fields of the lm class
#'
#' @importFrom methods new
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @export
slm.class <- setClass("slm",
         slots=list(method_cov_st="character",
                    cov_st = "numeric",
                    Cov_ST = "matrix",
                    model_selec = "numeric",
                    norm_matrix = "matrix",
                    design_qr = "matrix"),
         contains = "lm"
        )

#' @title Fitting Stationary Linear Models
#'
#' @description \code{slm} is used to fit linear models when the error process is assumed to be strictly stationary.
#'
#' @param myformula an object of class "\code{\link[stats:formula]{formula}}" (or one that can be coerced to that class):
#'  a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by \code{\link[base:as.data.frame]{as.data.frame}} to a data frame)
#'  containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)},
#'  typically the environment from which \code{slm} is called.
#' @param model,x,y,qr logicals. If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the response, the QR decomposition) are returned.
#' @param method_cov_st the method chosen by the user to estimate the autocovariances of the error process. The user
#'  has the choice between the methods "fitAR", "spectralproj", "efromovich", "kernel", "select" or "hac". By default, the
#'  "fitAR" method is used.
#' @param cov_st the estimated autocovariances of the error process. The user can give his own vector.
#' @param Cov_ST an argument given by the user if he wants to use his own covariance matrix estimator.
#' @param model_selec the order of the method. If \code{model_selec = -1}, the method works automatically.
#(depends on the method : AIC for AR, ...), model chosen by the user, -1 for automatic selection.
#' @param model_max maximal order of the method.
#' @param kernel_fonc to use if \code{method_cov_st = kernel}. Define the kernel to use in the method. The user can give his own kernel function.
#' @param block_size size of the bootstrap blocks if \code{method_cov_st = kernel}. \code{block_size} must be greater than \code{model_max}.
#' @param block_n blocks number to use for the bootstrap if \code{method_cov_st = kernel}.
#' @param plot logical. By default, \code{plot = FALSE}.
#'
#' @details The \code{slm} function is based on the architecture of the \code{lm} function.
#'  Models for \code{slm} are specified symbolically.
#'  A typical model has the form \code{response ~ terms} where \code{response} is the (numeric)
#'  response vector and \code{terms} is a series of terms which specifies a linear predictor for \code{response}.
#'  See the documentation of \code{\link[stats:lm]{lm}} for more details.
#'
#' @return \code{slm} returns an object of \code{\link[base:class]{class}} "\code{slm}". The function
#' \code{summary} is used to obtain and print a summary of the results. The generic accessor functions
#' \code{coefficients}, \code{effects}, \code{fitted.values} and \code{residuals} extract various
#'  useful features of the value returned by \code{slm}.
#'  An object of class "\code{slm}" is a list containing at least the following components:
#'  \item{method_cov_st }{ print the method chosen.}
#'  \item{cov_st }{ the estimated autocovariances of the error process. NA if "hac" is used.}
#'  \item{Cov_ST }{ if given by the user, the estimated covariance matrix of the error process. NA if "hac" is used.}
#If \code{method_cov_st = "hac"}, then \code{Cov_ST} is the covariance matrix estimator.}
#'  \item{model_selec }{ the order of the method.}
#'  \item{norm_matrix }{ the normalization matrix of the least squares estimator.}
#'  \item{design_qr }{ the matrix \eqn{(X^{t} X)^{-1}}.}
#'  \item{coefficients }{ a named vector of the estimated coefficients.}
#'  \item{residuals }{ the residuals, that is response minus fitted values.}
#'  \item{fitted.values }{ the fitted values.}
#'  \item{rank }{ the numeric rank of the fitted linear model.}
#'  \item{df.residual}{ the number of observations minus the number of variables.}
#'  \item{call }{ the matched call.}
#'  \item{terms }{ the \code{\link[stats:terms]{terms}} object used.}
#\item{contrasts }{ (only where relevant) the contrasts used.}
#'  \item{xlevels }{ (only where relevant) a record of the levels of the factors used in fitting.}
#\item{offset }{ the offset used (missing if none were used).}
#'  \item{y }{ if requested, the response used.}
#'  \item{x }{ if requested, the model matrix used.}
#'  \item{model }{ if requested (the default), the model frame used.}
#'
#' @export
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @seealso
#' \code{\link[slm:summary.slm]{summary}} for summaries.
#'
#' The generic functions \code{\link[stats:coef]{coef}}, \code{\link[stats:effects]{effects}},
#' \code{\link[stats:residuals]{residuals}}, \code{\link[stats:fitted]{fitted}}, \code{\link[stats:vcov]{vcov}}.
#'
#' \code{\link[slm:predict.slm]{predict}} for prediction, including confidence intervals for \eqn{x' beta}, where \eqn{x'} is a new observation of the design.
#'
#' \code{\link[slm:confint.slm]{confint}} for confidence intervals of \emph{parameters}.
#'
#' @examples
#' data("shan")
#' slm(shan$PM_Xuhui ~ . , data = shan, method_cov_st = "fitAR", model_selec = -1)
#'
#' data("co2")
#' y = as.vector(co2)
#' x = as.vector(time(co2)) - 1958
#' reg1 = slm(y ~ x + I(x^2) + I(x^3) + sin(2*pi*x) + cos(2*pi*x) + sin(4*pi*x) +
#'  cos(4*pi*x) + sin(6*pi*x) + cos(6*pi*x) + sin(8*pi*x) + cos(8*pi*x),
#'  method_cov_st = "fitAR", model_selec = -1, plot = TRUE)
#'
#' reg2 = slm(y ~ x + I(x^2) + I(x^3) + sin(2*pi*x) + cos(2*pi*x) + sin(4*pi*x) +
#'  cos(4*pi*x) + sin(6*pi*x) + cos(6*pi*x) + sin(8*pi*x) + cos(8*pi*x),
#'  method_cov_st = "kernel", model_selec = -1, model_max = 50, kernel_fonc = triangle,
#'  block_size = 100, block_n = 100)
slm <- function(myformula,
                data = NULL,
                model = TRUE,
                x = FALSE,
                y = FALSE,
                qr = TRUE,
                method_cov_st="fitAR",
                cov_st = NULL,
                Cov_ST = NULL,
                model_selec = -1,
                model_max = 50,
                kernel_fonc = NULL,
                block_size = NULL,
                block_n = NULL,
                plot = FALSE){

  #least square estimation from lm procedure
  lm_call <- lm(myformula, data = data, model = model, x = x, y = y, qr = qr)
  lm_call$call = "slm(formula = myformula, data = data, x = x, y = y)"

  #normalization matrix
  Y = lm_call$model[[1]]
  #intercept added in the design X. (To compute Cn)
  design = model.matrix(lm_call)

  p <- lm_call$rank #p = number of variables, with intercept
  Qr <- qr(lm_call)
  p1 <- 1L:p #p = var. number
  design_qr <- chol2inv(Qr$qr[p1, p1, drop = FALSE]) # = (X^{t} X)^{-1}, with intercept
  #design_qr = chol2inv(qr.R(qr(design)))

  #norm_matrix = diag(sqrt(diag(solve(vcov(lm_call)/summary(lm_call)$sigma^2))))
  norm_matrix = diag(sqrt(apply(design^2,2,sum)), nrow=dim(design)[2])

  # covariance estimation
  if (is.null(cov_st) && is.null(Cov_ST)){
    if (method_cov_st=="hac") {
      model_selec = NA_real_
      cov_st = NA_real_
      Cov_ST = matrix(NA_real_)
    } else {
      mylist = cov_method(epsilon = lm_call$residuals,
                          method_cov_st = method_cov_st,
                          model_selec = model_selec,
                          model_max = model_max,
                          kernel_fonc = kernel_fonc,
                          block_size = block_size,
                          block_n = block_n,
                          plot = plot)
      cov_st = mylist$cov_st
      Cov_ST = matrix(NA_real_)
      model_selec = mylist$model_selec
    }
  } else {
    if (is.null(Cov_ST)) {
      epsilon = lm_call$residuals
      method_cov_st = "manual"
      model_selec = NA_real_
      cov_st = cov_st
      Cov_ST = matrix(NA_real_)
    } else if (is.null(cov_st)) {
      epsilon = lm_call$residuals
      method_cov_st = "manual_matrix"
      model_selec = NA_real_
      cov_st = NA_real_
      Cov_ST = Cov_ST
    } else {
      epsilon = lm_call$residuals
      method_cov_st = "manual_matrix"
      model_selec = NA_real_
      cov_st = NA_real_
      Cov_ST = Cov_ST #if both, keep matrix
    }
  }

  # return an object of class slm
  out <- slm.class(lm_call,
                   method_cov_st=method_cov_st,
                   cov_st = cov_st,
                   Cov_ST = Cov_ST,
                   model_selec = model_selec,
                   norm_matrix = norm_matrix,
                   design_qr = design_qr )

  return(out)
}
#cov_st = covariance vector of the error process (hat_epsilon, residual based)
#Cov_ST = covariance matrix of the error process (hat_espilon, residual based)

#' @title Methods to estimate the autocovariances of a process
#'
#' @description This function gives the estimation of the autocovariances of the error process, with the method chosen by the user.
#'  Five methods are available: "fitAR", "spectralproj", "efromovich", "kernel" and "select".
#'
#' @param method_cov_st the method chosen by the user.
#' @param epsilon an univariate process.
#' @param model_selec the order of the method. If \code{model_selec = -1}, the method works automatically.
#(depends on the method : AIC for AR, ...)
#' @param model_max maximal dimension of the method.
#' @param kernel_fonc to use if \code{method_cov_st = kernel}. Define the kernel to use in the method. The user can give his own kernel function.
#' @param block_size size of the bootstrap blocks if \code{method_cov_st = kernel}. \code{block_size} must be greater than \code{model_max}.
#' @param block_n blocks number to use for the bootstrap if \code{method_cov_st = kernel}.
#' @param plot logical. By default, \code{plot = FALSE}.
#'
#' @return The function returns the autocovariances computed with the chosen method.
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @export
#'
#' @examples
#' x = arima.sim(list(ar=c(0.4,0.2)),1000)
#' cov_method(x, method_cov_st = "fitAR", model_selec = -1)
cov_method <- function(epsilon, method_cov_st = "fitAR", model_selec = -1,
                       model_max=NULL, kernel_fonc = NULL, block_size = NULL,
                       block_n = NULL, plot = FALSE){

  switch(method_cov_st,
         # fitting an AR on the residuals
         fitAR={out = cov_AR(epsilon,model_selec = model_selec,plot=plot)},
         # via spectral density estimation
         spectralproj={out = cov_spectralproj(epsilon,model_selec = model_selec,model_max = model_max,plot=plot)},
         # efromovich method
         efromovich={out = cov_efromovich(epsilon,plot=plot)},
         # kernel method
         kernel={out = cov_kernel(epsilon,model_selec = model_selec,model_max = model_max,kernel_fonc = kernel_fonc
                                  ,block_size = block_size,block_n = block_n,plot=plot)},
         # manual
         select={out = cov_select(epsilon,model_selec = model_selec,plot=plot)}
         )

  return(out)
}

#' @title Covariance estimation by AR fitting
#'
#' @description Fit an autoregressive model to the process and compute the theoretical autocovariances of the fitted AR process.
#'  By default, the order is chosen by using the AIC criterion (\code{model_selec = -1}).
#'
#' @param epsilon an univariate process.
#' @param model_selec the order of the method. If \code{model_selec = -1}, it is chosen automatically by using the AIC criterion.
#' @param plot logical. By default, \code{plot = FALSE}. If \code{plot = TRUE}, then the ACF and the PACF of the vector \code{epsilon} is plotted.
#'
#' @return The function returns the vector of the theoretical autocovariances of the AR process fitted on the process \code{epsilon}.
#'  \item{model_selec}{the order selected.}
#'  \item{cov_st}{the vector of theoretical autocovariances of the fitted AR process.}
#'
#' @importFrom stats acf pacf ar
#'
#' @references
#'  P.J. Brockwell and R.A. Davis (1991). Time Series: Theory and Methods. \emph{Springer Science & Business Media}.
#'
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @export
#'
#' @examples
#' x = arima.sim(list(ar=c(0.4,0.2)),1000)
#' cov_AR(x, model_selec = 2, plot = TRUE)
cov_AR <- function(epsilon, model_selec = -1, plot=FALSE){
  if (plot == TRUE) {
    acf(epsilon, lag.max=sqrt(length(epsilon)), type="correlation", main=" ")#ACF of the residuals")
    pacf(epsilon, lag.max=sqrt(length(epsilon)), main=" ")#PACF of the residuals")
  }
  n = length(epsilon)
  if (model_selec == -1){
    my_ar = ar(epsilon, aic = TRUE)
    if (my_ar$order==0) {
      cov_st = rep(0,n)
      cov_st[1] = var(epsilon)
      model_selec = my_ar$order #order of the ar process
    } else {
      coef_ar = my_ar$ar #coefficients of the ar process
      cov_st = ltsa::tacvfARMA(phi=coef_ar, maxLag=n-1, sigma2=my_ar$var.pred) #theoretical covariance vector of the ar fitted
      model_selec = my_ar$order #order of the ar process
    }
  } else {
    if (model_selec == 0){
      cov_st = rep(0,n)
      cov_st[1] = var(epsilon)
    } else {
      my_ar = ar(epsilon, aic = FALSE, order.max = model_selec)
      coef_ar = my_ar$ar #coefficients of the ar process
      cov_st = ltsa::tacvfARMA(phi=coef_ar, maxLag=n-1, sigma2=my_ar$var.pred) #theoretical covariance vector of the ar fitted
    }
  }
  return(list(model_selec=model_selec,cov_st=cov_st))
}

#' @title Data-driven spectral density estimation
#'
#' @description Computes a data-driven histogram estimator of the spectral density of a process and compute its Fourier coefficients,
#'  that is the associated autocovariances. For a dimension \eqn{d}, the estimator of the spectral density is an histogram on a regular basis of
#'  size \eqn{d}. Then we use a penalized criterion in order to choose the dimension which balance the bias and the variance, as proposed in Comte (2001). The penalty
#'  is of the form \eqn{c*d/n}, where \eqn{c} is the constant and \eqn{n} the sample size. The dimension and the constant of the penalty are
#'  chosen with the slope heuristic method, with the dimension jump algorithm (from package "\code{\link[capushe:capushe]{capushe}}").
#'
#' @usage cov_spectralproj(epsilon, model_selec = -1,
#'  model_max = min(100,length(epsilon)/2), plot = FALSE)
#'
#' @param epsilon an univariate process.
#' @param model_selec the dimension of the method. If \code{model_selec = -1}, the method works automatically and take a dimension between 1 and \code{model_max}.
#' @param model_max the maximal dimension. By default, it is equal to the minimum between 100 and the length of the process divided by 2.
#' @param plot logical. By default, \code{plot = FALSE}. If \code{plot = TRUE}, plot the spectral density estimator of the process.
#'
#' @return The function returns the estimated autocovariances of the process, that is the Fourier coefficients
#'  of the spectral density estimates, and the dimension chosen by the algorithm.
#'  \item{model_selec}{the dimension selected.}
#'  \item{cov_st}{the estimated autocovariances.}
#'
#' @export
#'
#' @references
#'  J.P. Baudry, C. Maugis B. and Michel (2012). Slope heuristics: overview and implementation. \emph{Statistics and Computing}, 22(2), 455–470.
#'
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#'  F. Comte (2001). Adaptive estimation of the spectrum of a stationary Gaussian sequence. \emph{Bernoulli}, 7(2), 267-298.
#'
#' @seealso The R package \code{\link[capushe]{capushe}}.
#'
#' @seealso Slope heuristic algorithm \code{\link[capushe]{DDSE}}.
#'
#' @seealso Dimension jump algorithm \code{\link[capushe]{Djump}}.
#'
#' @examples
#' x = arima.sim(list(ar=c(0.2), ma=c(0.3,0.05)), n=100)
#' cov_spectralproj(x, model_selec = -1)
cov_spectralproj <- function(epsilon, model_selec = -1, model_max = min(100,length(epsilon)/2), plot=FALSE){
  n = length(epsilon)
  cov_epsilon = as.vector(acf(epsilon,type="covariance",lag.max = n-1,plot=FALSE)$acf) #covariance vector of the residuals

  if (model_selec==-1) {
    mat_a = matrix(0, nrow=model_max, ncol=model_max) #matrix with the estimated coefficients of the spectral density, for every dimension d
    contrast = rep(0,model_max) #vector of contrasts
    pen = rep(0,model_max) #vector of penalties
    pen_contrast = rep(0,model_max) #vector of the penalized contrasts

    for(d in seq(1,model_max)) {
      a_hat = rep(0,d) #vector of the estimated coefficients, for a fixed dimension d
      for (j in seq(1,d)) {
        #j = 0, ..., d-1
        vec = rep(0,n-1)
        for (r in seq(1,n-1)) {
          vec[r] = (cov_epsilon[r+1]/r)*(sin((pi*j*r)/d) - sin((pi*(j-1)*r)/d)) #for the sum in the expression of the estimated coefficients
        }
        a_hat[j] = sqrt(d/pi)*(cov_epsilon[1]/(2*d) + (1/pi)*sum(vec))
      }
      mat_a[d,1:d] = a_hat #for a fixed dimension d, raw vector of the estimated coefficients
      pen[d] = d #penalty, just Dm=d
      contrast[d] = (-1)*sum(mat_a[d,1:d]^2)
      pen_contrast[d] = contrast[d] + pen[d] #penalized contrast without the constant
    }

    # dimension selection
    datacap = matrix(0,nrow=model_max,ncol=4,dimnames = list(seq(1,model_max),c("model","pen","complexity","contrast")))
    datacap[,1] = pen
    datacap[,2] = pen
    datacap[,3] = pen
    datacap[,4] = contrast

    #with slope heuristic
    #d_hat = capushe::DDSE(datacap)

    #with dimension jump
    d_hat = capushe::Djump(datacap)

    dhat = as.numeric(d_hat@model)
    model_selec = dhat
    spec_dens = sqrt(dhat/pi)*mat_a[dhat,(1:dhat)] #spectral density estimator in [0,pi], for the chosen dimension
  } else {
    dhat = model_selec
    a_hat = rep(0,dhat)
    for (j in seq(1,dhat)) {
      #j = 0, ..., d-1
      vec = rep(0,n-1)
      for (r in seq(1,n-1)) {
        vec[r] = (cov_epsilon[r+1]/r)*(sin((pi*j*r)/dhat) - sin((pi*(j-1)*r)/dhat))
      }
      a_hat[j] = sqrt(dhat/pi)*(cov_epsilon[1]/(2*dhat) + (1/pi)*sum(vec))
    }
    spec_dens = sqrt(dhat/pi)*a_hat #spectral density estimator in [0,pi], for the chosen dimension
  }

  if (plot==TRUE) {
    x = seq(0,pi-(pi/dhat),by=pi/dhat) #def of the histogram basis
    plot(x,spec_dens,type="s",ylab="spectral density")
  }

  #Fourier coefficients
  cov_st = rep(0,n)
  cov_st[1] = ((2*pi)/dhat)*sum(spec_dens)
  interm = rep(0,dhat)
  for (k in seq(2,n)) {
    for (j in seq(1,dhat)) {
      interm[j] = spec_dens[j]*(sin(((k-1)*pi*j)/dhat) - sin(((k-1)*pi*(j-1))/dhat))
    }
    cov_st[k] = (2/(k-1))*sum(interm)
  }
  return(list(model_selec=model_selec,cov_st=cov_st))
}

# 0 -> variance
# 1 -> gamma(1), etc.
#' @title Covariances Selection
#'
#' @description Allows the user to select the lags of the autocovariance terms of the process to be kept.
#'
#' @param epsilon an univariate process.
#' @param model_selec a vector with the positive lags of the selected autocovariance terms. The variance (lag = 0) is automatically selected.
#' @param plot logical. By default, \code{plot = FALSE}. If \code{plot = TRUE} the ACF of the process is plotted.
#'
#' @details In the framework of \code{slm}, this is a manual method for estimating the covariance matrix of the error process
#'  by only selecting some autocovariance terms from the residual autocovariances.
#'
#' @return This function returns the estimated autocovariance terms.
#'  \item{model_selec}{the vector with the positive lag of the selected autocovariance terms.}
#'  \item{cov_st}{the vector of the selected autocovariances.}
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @export
#'
#' @examples
#' x = arima.sim(list(ar=c(0.2,0.1,0.25)),1000)
#' cov_select(x, c(1,3,5))
cov_select <- function(epsilon, model_selec, plot=FALSE){
  n = length(epsilon)
  if (plot==TRUE) {
    acf(epsilon,type="correlation",lag.max=max(model_selec)+1,main=" ")#"ACF of the residuals")
  }
  cov_epsilon = as.vector(acf(epsilon,type="covariance",lag.max = n-1,plot=FALSE)$acf)
  cov_st = rep(0,n)
  cov_st[1] = cov_epsilon[1]
  cov_st[model_selec+1] = cov_epsilon[model_selec+1]
  return(list(model_selec=model_selec,cov_st=cov_st))
}

#' @title Kernel estimation: bootstrap method
#'
#' @description This method estimates the spectral density and the autocovariances of the error process via a lag-window
#'  (or kernel) estimator (see P.J. Brockwell and R.A. Davis (1991). Time Series: Theory and Methods. \emph{Springer Science & Business Media},
#'  page 330). The weights are computed according to a kernel \eqn{K} and a bandwidth \eqn{h} (or a lag),
#'  to be chosen by the user. The lag can be computed automatically by using a bootstrap technique (as in Wu and Pourahmadi (2009)), via the \code{\link[slm:Rboot]{Rboot}} function.
#which can be understood like the number of covariance terms to keep to have a good approximation of the covariance vector.
#'
#' @usage cov_kernel(epsilon, model_selec = -1,
#'  model_max = min(50,length(epsilon)/2), kernel_fonc = triangle,
#'  block_size = length(epsilon)/2, block_n = 100, plot = FALSE)
#'
#' @param epsilon an univariate process.
#' @param model_selec the order of the method. If \code{model_selec = -1}, the method chooses
#'  the treshold automatically. If \code{model_selec = k}, then only \code{k} autocovariance terms are kept
#'  and smoothed by the kernel.
#' @param model_max the maximal order.
#' @param kernel_fonc define the kernel to use in the method. The user can give his own kernel function.
#' @param block_size size of the bootstrap blocks. \code{block_size} must be greater than \code{model_max}.
#' @param block_n blocks number to use for the bootstrap.
#' @param plot logical. By default, \code{plot = FALSE}. If \code{plot = TRUE}, the risk curve is returned and the
#'  ACF of the process.
#'
#' @return The method returns the tapered autocovariance vector with \code{model_selec} autocovariance terms.
#'  \item{model_selec}{the number of selected autocovariance terms.}
#'  \item{cov_st}{the estimated autocovariances.}
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#'  W.B. Wu, M. Pourahmadi (2009). Banding sample autocovariance matrices of stationary processes. \emph{Statistica Sinica}, pp. 1755–1768.
#'
#' @export
#'
#' @examples
#' x = arima.sim(list(ar=c(0.7)),1000)
#' cov_kernel(x, model_selec = -1, block_n = 10, plot = TRUE)
cov_kernel <- function(epsilon, model_selec = -1, model_max = min(50,length(epsilon)/2), kernel_fonc = triangle,
                       block_size = length(epsilon)/2, block_n = 100, plot=FALSE){
  n = length(epsilon)
  if (model_selec==-1) {
    #bootstrap, one se rule
    risk = rep(0,model_max) #risk
    SE = rep(0,model_max) #bootstrap standard error for every lag
    #treshold from 1 (only gamma0) to model_max (from gamma0 to gamma(model_max-1))
    for (treshold in seq(1,model_max)) {
      result = Rboot(epsilon,treshold,block_size,block_n,model_max,kernel_fonc)
      risk[treshold] = result[1]
      SE[treshold] = result[2]
      #risk curve with bootstrap method
    }
    model_selec = which.min(risk)
    if (plot==TRUE) {
      plot(seq(0,model_max-1),risk,type="l",xlab="lag")
      acf(epsilon,type="correlation",lag.max=model_selec-1,main=" ")
    }
    cov_st = rep(0,n)
    cov_st[1:model_selec] = acf(epsilon,type="covariance",lag.max=model_selec-1,plot=FALSE)$acf
    kern = rep(0,n)
    kern[1:model_selec] = kernel_fonc((0:(model_selec-1))/model_selec)
    cov_st = kern*cov_st
    return(list(model_selec=model_selec-1,cov_st=cov_st))
  } else {
    if (plot==TRUE) {
      acf(epsilon,lag.max=sqrt(length(epsilon)),type="correlation", main=" ")#ACF of the residuals")
    }
    model_selec = model_selec + 1
    cov_st = rep(0,n)
    cov_st[1:model_selec] = acf(epsilon,type="covariance",lag.max=model_selec-1,plot=FALSE)$acf
    kern = rep(0,n)
    kern[1:model_selec] = kernel_fonc((0:(model_selec-1))/model_selec)
    cov_st = kern*cov_st
    return(list(model_selec=model_selec-1,cov_st=cov_st))
  }
}

#model_selec=-1 auto
#model_selec=0 keep only variance
#model_selec=k keep variance + covariances up to lag k

#' @title Spectral density estimation: Efromovich method
#'
#' @description This method estimates the spectral density and the autocovariances of the error process via a lag-window
#'  estimator based on the rectangular kernel (see P.J. Brockwell and R.A. Davis (1991). Time Series: Theory and Methods.
#'  \emph{Springer Science & Business Media}, page 330). The lag is computed according to Efromovich's algorithm (Efromovich (1998)).
#and we use its Fourier coefficients to build an estimation of the covariance vector of the process.
#'
#' @param epsilon an univariate process.
#' @param plot logical. By default, \code{plot = FALSE}. If \code{plot = TRUE}, the ACF of the process \code{epsilon}
#'  is plotted.
#'
#' @return The function returns the estimated autocovariances of the process, that is the Fourier coefficients
#'  of the spectral density estimates, and the order chosen by the algorithm.
#'  \item{model_selec}{the number of selected autocovariance terms.}
#'  \item{cov_st}{the estimated autocovariances.}
#'
#' @export
#'
#' @references
#'  P.J. Brockwell and R.A. Davis (1991). Time Series: Theory and Methods. \emph{Springer Science & Business Media}.
#'
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#'  S. Efromovich (1998). Data-driven efficient estimation of the spectral density. \emph{Journal of the American Statistical Association}, 93(442), 762-769.
#'
#' @examples
#' x = arima.sim(list(ar=c(0.4,0.2)),1000)
#' cov_efromovich(x)
#algo small samples for n<=10000
cov_efromovich <- function(epsilon, plot = FALSE) {
  n = length(epsilon)
  #Step 1
  Jn = floor((log(n))^(5/4))
  cov_epsilon = as.vector(acf(epsilon,type="covariance",lag.max=n-1,plot=FALSE)$acf) #covariance vector of the residuals
  #Step 2
  dn_hat = cov_epsilon[1]^2 + 2*sum(cov_epsilon[2:(Jn+1)]^2)
  #Step 3
  Gamma_hat = rep(0,(Jn+1))
  for (j in seq(1:(Jn+1))) {
    Gamma_hat[j] = max(0,(cov_epsilon[j]^2 - dn_hat/n))
  }
  wn_hat = rep(0,n)
  #wn_hat = rep(0,(Jn+1))
  for (j in seq(1:(Jn+1))) {
    wn_hat[j] = Gamma_hat[j]/(cov_epsilon[j]^2)
  }
  #Step 4
  sum_jn = rep(0,(Jn+1))
  sum_jn[1] = 2*dn_hat/n - (cov_epsilon[1])^2 #case J=0
  for (j in seq(2,(Jn+1))) {
    sum_jn[j] = sum_jn[j-1] + 2*(2*dn_hat/n - cov_epsilon[j]^2)
  }
  Jn_hat = which.min(sum_jn)

  wn_hat_2 = rep(0,n)
  for (j in seq(1:(Jn_hat))) {
    wn_hat_2[j] = wn_hat[j]
  }
  #Step 5:
  cov_st = wn_hat_2*cov_epsilon
  model_selec = Jn_hat - 1
  if (plot==TRUE) {
    acf(epsilon,type="correlation",lag.max=model_selec,main=" ")#ACF of the residuals")
  }
  return(list(model_selec=model_selec,cov_st=cov_st))
  # } else {
  #   bn = 1/(log(log(n+20)))
  #   cov_epsilon = as.vector(acf(epsilon,type="covariance",lag.max=n-1)$acf) #covariance vector of the residuals
  #   b_inf = floor(bn*log(n))
  #   b_sup = floor(log(n)/bn)
  #   dn_hat = cov_epsilon[1]^{2} + 2*sum(cov_epsilon[2:(floor(b_sup)+1)]^2)
  #   B_hat = rep(0,(b_sup-b_inf+1))
  #   J_hat = 0
  #   for (J in seq(b_inf,b_sup)) {
  #     B_hat[J-b_inf+1] = (1/pi)*(sum(cov_epsilon[(1+J):(b_sup)]^2) - (b_sup-J)*(dn_hat/n))
  #     if (B_hat[J-b_inf+1] < ((bn*log(n))/n)) {
  #       J_hat=J
  #       break
  #     }
  #   }
  #   cov_epsilon[(J_hat+2):n] = 0
  #   cov_st = cov_epsilon
  #   model_selec = J_hat
  #   if (plot==TRUE) {
  #     acf(epsilon,type="correlation",lag.max=model_selec,main=" ")#ACF of the residuals")
  #   }
  #   return(list(model_selec=model_selec,cov_st=cov_st))
  # }
}


# F test functions ?

# plot functions

# predict functions

# un test anova specifique ?

# backward-stepwise functions ?


