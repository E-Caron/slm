#' @title Confidence intervals for the Model Parameters
#'
#' @description Computes confidence intervals for the model parameters.
#'
#' @param object a fitted model object of class \code{slm}.
#' @param parm 	a specification of which parameters are to be given confidence intervals,
#'  that is a vector of numbers. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... additional argument(s) for methods.
#'
#' @return This function returns the confidence intervals for the parameters of the model.
#'
#' @importFrom stats lm confint qnorm
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @seealso
#' \code{\link[stats:confint.lm]{confint.lm}}.
#'
#' @export
#'
#' @examples
#' data("shan")
#' reg1 = slm(shan$PM_Xuhui ~ . , data = shan, method_cov_st = "fitAR", model_selec = -1)
#' confint(reg1, level = 0.8)
#'
#' data("co2")
#' y = as.vector(co2)
#' x = as.vector(time(co2)) - 1958
#' reg2 = slm(y ~ x + I(x^2) + I(x^3) + sin(2*pi*x) + cos(2*pi*x) + sin(4*pi*x) +
#'  cos(4*pi*x) + sin(6*pi*x) + cos(6*pi*x) + sin(8*pi*x) + cos(8*pi*x),
#'  method_cov_st = "fitAR", model_selec = -1, plot = TRUE)
#' confint(reg2, level = 0.9)
confint.slm = function(object, parm=NULL, level=0.95, ...) {
  mod = object$model
  names(mod)[1] = "Y"
  reglm = lm(Y~., data=mod)
  norm_matrix = diag(object@norm_matrix)
  coefs = object$coefficients
  Cn = diag(cov_matrix_estimator(object))
  levels = 1-level
  if (is.null(parm)) {
    interval = confint(reglm, level=level, ...)
    interval[1:dim(interval)[1],1] = coefs - (qnorm(1-levels/2)*sqrt(Cn))/norm_matrix
    interval[1:dim(interval)[1],2] = coefs + (qnorm(1-levels/2)*sqrt(Cn))/norm_matrix
  } else {
    interval = confint(reglm, parm=parm, level=level, ...)
    interval[1:dim(interval)[1],1] = coefs[parm] - (qnorm(1-levels/2)*sqrt(Cn[parm]))/norm_matrix[parm]
    interval[1:dim(interval)[1],2] = coefs[parm] + (qnorm(1-levels/2)*sqrt(Cn[parm]))/norm_matrix[parm]
  }
  return(interval)
}
# confint.slm = function(object, parm=NULL, level=0.95, ...) {
#   Y = as.matrix(object$model[1])
#   design = cbind(rep(1,length(Y)),as.matrix(object$model[-1])) #intercept added in the design X.
#   num_var = dim(design)[2]
#   norm_matrix = diag(object@norm_matrix)
#   coefs = object$coefficients
#   Cn = diag(cov_matrix_estimator(object))
#   levels = 1-level
#   if (is.null(parm)) {
#     interval = matrix(0, nrow=num_var, ncol=2, dimnames = list(c("(Intercept)",names(object$model[,-1])),c(levels/2, 1-levels/2)))
#     interval[1:num_var,1] = coefs - (qnorm(1-levels/2)*sqrt(Cn))/norm_matrix
#     interval[1:num_var,2] = coefs + (qnorm(1-levels/2)*sqrt(Cn))/norm_matrix
#   } else {
#     interval = matrix(0, nrow=length(parm), ncol=2, dimnames = list(parm,c(levels/2, 1-levels/2)))
#     interval[1:length(parm),1] = coefs[parm] - (qnorm(1-levels/2)*sqrt(Cn[parm]))/norm_matrix[parm]
#     interval[1:length(parm),2] = coefs[parm] + (qnorm(1-levels/2)*sqrt(Cn[parm]))/norm_matrix[parm]
#   }
#   return(interval)
# }

#' @title Summarizing Stationary Linear Model Fits
#'
#' @description Summary method for class "\code{slm}".
#'
#' @param object an object of class "\code{slm}", usually, a result of a call to \code{slm}.
#' @param correlation logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed.
#' @param symbolic.cor logical. If \code{TRUE}, print the correlations in a symbolic form (see \code{\link[stats:symnum]{symnum}}) rather than as numbers.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The function \code{summary.slm} computes and returns a list of summary statistics of the fitted linear model
#'  given in \code{object}, using the components (list elements) "\code{call}" and "\code{terms}" from its argument, plus:
#'  \item{residuals }{ the residuals, that is response minus fitted values.}
#'  \item{coefficients }{ a \eqn{p*4} matrix with columns for the estimated coefficient, its standard error, z-statistic and corresponding (two-sided) p-value.
#'   Aliased coefficients are omitted.}
#'  \item{aliased }{ named logical vector showing if the original coefficients are aliased.}
#'  \item{sigma }{ the square root of the estimated variance of the error process.}
#'  \item{df }{ degrees of freedom, a 3-vector \eqn{(p, n-p, p*)}, the first being the number of non-aliased coefficients, the last being the total number of coefficients.}
#'  \item{chi2statistic }{ a 2-vector with the value of the chi2-statistic with its degree of freedom.}
#'  \item{r.squared }{ \eqn{R^2}, the 'fraction of variance explained by the model'.}
#'  \item{cov.unscaled }{ a \eqn{p*p} matrix of (unscaled) covariances of the \eqn{coef[j], j=1, \ldots, p}.}
#'  \item{correlation }{ the correlation matrix corresponding to the above \code{cov.unscaled}, if \code{correlation = TRUE} is specified.}
#'  \item{symbolic.cor }{ (only if \code{correlation} is true.) The value of the argument \code{symbolic.cor}.}
#'
#' @import stats
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @seealso
#'  The model fitting function \code{\link[slm:slm]{slm}}, \code{\link[base:summary]{summary}}.
#'
#'  The function \code{\link[stats:coef]{coef}} extracts the matrix of coefficients with standard errors, z-statistics and p-values.
#'
#' @export
#'
#' @examples
#' data("shan")
#' reg1 = slm(shan$PM_Xuhui ~ . , data = shan, method_cov_st = "fitAR", model_selec = -1)
#' summary(reg1)
#'
#' data("co2")
#' y = as.vector(co2)
#' x = as.vector(time(co2)) - 1958
#' reg2 = slm(y ~ x + I(x^2) + I(x^3) + sin(2*pi*x) + cos(2*pi*x) + sin(4*pi*x) +
#'  cos(4*pi*x) + sin(6*pi*x) + cos(6*pi*x) + sin(8*pi*x) + cos(8*pi*x),
#'  method_cov_st = "fitAR", model_selec = -1)
#' summary(reg2)
summary.slm <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  z <- object
  p <- z$rank #p = number of variables, with intercept
  rdf <- z$df.residual #n-p
  #case p==0
  if (p == 0) {
    r <- z$residuals #hat espilon
    n <- length(r) #obs. number
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    } else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf #resvar = sigma2 = rss/n-p = variance
    ans <- z[c("call", "terms", if(!is.null(z$weights)) "weights")]
    class(ans) <- "summary.slm"
    ans$aliased <- is.na(coef(object))  # used in print method
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <-
      list(NULL, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    ans$sigma <- sqrt(resvar) #Residual Standard Error
    ans$r.squared <- 0
    #  ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms))
    stop("invalid 'slm' object:  no 'terms' component")
  if(!inherits(object, "slm"))
    warning("calling summary.slm(<fake-slm-object>) ...")
  Qr <- qr.slm(object)
  n <- NROW(Qr$qr) #obs. number
  if(is.na(z$df.residual) || n - p != z$df.residual)
    warning("residual degrees of freedom in object suggest this is not an \"slm\" fit")
  ## do not want missing values substituted here
  r <- z$residuals #residuals hat epsilon
  f <- z$fitted.values #hat Y
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept"))
      sum((f - mean(f))^2) else sum(f^2)
    rss <- sum(r^2)
  } else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f /sum(w))
      sum(w * (f - m)^2)
    } else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  resvar <- rss/rdf #sigma2
  ## see thread at https://stat.ethz.ch/pipermail/r-help/2014-March/367585.html
  if (is.finite(resvar) &&
      resvar < (mean(f)^2 + var(f)) * 1e-30)  # a few times .Machine$double.eps^2
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- 1L:p #p = var. number

  #change table values
  norm_matrix = object@norm_matrix
  num_var = dim(norm_matrix)[2]

  #Cn: covariance matrix estimator
  Cn = cov_matrix_estimator(object)

  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE]) #= (XtX)^-1

  #Modification of the estimation of the standard error of each beta_chapeau_j
  se <- sqrt(diag(Cn))/diag(norm_matrix) #sqrt(diag(R) * resvar) #Std. Error of hat beta j
  est <- z$coefficients[Qr$pivot[p1]] #Estimate
  #Modification of the t-value (-> z-value). Student's test H0 : beta_j = 0.
  zval <- est/se #tval <- est/se #"new Student" tests
  ans <- z[c("call", "terms", if(!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <-
    #cbind(est, se, zval, 2*pnorm(abs(zval), mean=0, sd=1, lower.tail=FALSE))
    #Modification of the p-value, stars are automatically modified
    cbind(est, se, zval, 2*(1-pnorm(abs(zval), mean=0, sd=1))) #2*pt(abs(tval), rdf, lower.tail = FALSE)) #p-value
  #attention aux p-values, apparemment ici il faut ecrire 1-pnorm... alors que pour la fstatistic, juste pchisq... et R fait
  #le reste automatiquement...?
  dimnames(ans$coefficients) <-
    list(names(z$coefficients)[Qr$pivot[p1]],
         c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  ans$aliased <- is.na(coef(object))  # used in print method
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf, NCOL(Qr$qr)) #p, n-p, p
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 1L else 0L
    ans$r.squared <- mss/(mss + rss)
    #ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
    #F-stat: full model vs intercept
    Z = solve(chol(Cn[2:num_var,2:num_var]))%*%(norm_matrix[2:num_var,2:num_var]%*%est[-1])
    ans$chi2statistic <- c(value = sum(Z^2), numdf = num_var - 1)

    #ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
    #                    numdf = p - df.int, dendf = rdf)
  } else ans$r.squared <- 0 #else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R #cov.unscaled = (XtX)^-1
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se) #outer = produit tensoriel
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if(!is.null(z$na.action)) ans$na.action <- z$na.action
  class(ans) <- "summary.slm"
  ans
}

#intercept added in the design X. To compute Cn

# #' @title Print summary for slm
# #'
# #' @param x an object of class "\code{summary.slm}", usually, a result of a call to \code{summary.slm}.
# #' @param digits the number of significant digits to use when printing.
# #' @param symbolic.cor logical. If TRUE, print the correlations in a symbolic form (see \code{\link[stats:symnum]{symnum}}) rather than as numbers.
# #' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
# #' @param ... further arguments passed to or from other methods.
# #'
# #' @import stats
# #'
# #' @seealso
# #' \code{\link[stats:summary.lm]{summary.lm}}.
# #'
#' @export
print.summary.slm <-
  function (x, digits = max(3L, getOption("digits") - 3L),
            symbolic.cor = x$symbolic.cor,
            signif.stars = getOption("show.signif.stars"),	...)
  {
    cat("\nCall:\n", # S has ' ' instead of '\n'
        paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep = "") #call display
    resid <- x$residuals #hat espilon
    df <- x$df #n-p
    rdf <- df[2L]
    cat(if(!is.null(x$weights) && diff(range(x$weights))) "Weighted ",
        "Residuals:\n", sep = "")
    if (rdf > 5L) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if (length(dim(resid)) == 2L)
        structure(apply(t(resid), 1L, quantile),
                  dimnames = list(nam, dimnames(resid)[[2L]]))
      else  {
        zz <- zapsmall(quantile(resid), digits + 1L)
        structure(zz, names = nam)
      }
      print(rq, digits = digits, ...)
    }
    else if (rdf > 0L) {
      print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
      cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
      cat("\n")
    }
    if (length(x$aliased) == 0L) {
      cat("\nNo Coefficients\n")
    } else {
      if (nsingular <- df[3L] - df[1L])
        cat("\nCoefficients: (", nsingular,
            " not defined because of singularities)\n", sep = "")
      else cat("\nCoefficients:\n")
      coefs <- x$coefficients
      if(!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4, dimnames=list(cn, colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
      }

      printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
    }
    ##
    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)))
    #    format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")
    if (!is.null(x$chi2statistic)) {
      cat("\nMultiple R-squared: ", formatC(x$r.squared, digits = digits))
      cat("\nchi2-statistic:", formatC(x$chi2statistic[1L], digits = digits), "on", x$chi2statistic[2L],"DF,  p-value:",
          format.pval(pchisq(x$chi2statistic[1L], x$chi2statistic[2L], lower.tail = FALSE), digits = digits))

      #cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
      #cat(",\tAdjusted R-squared: ",formatC(x$adj.r.squared, digits = digits),
      #    "\nF-statistic:", formatC(x$fstatistic[1L], digits = digits),
      #    "on", x$fstatistic[2L], "and",
      #    x$fstatistic[3L], "DF,  p-value:",
      #    format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
      #                   x$fstatistic[3L], lower.tail = FALSE),
      #                digits = digits))

      #cat("\tp-value:", formatC(x$fstatistic[1L], digits = digits))
      cat("\n")
    }
    correl <- x$correlation
    if (!is.null(correl)) {
      p <- NCOL(correl)
      if (p > 1L) {
        cat("\nCorrelation of Coefficients:\n")
        if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
          print(symnum(correl, abbr.colnames = NULL))
        } else {
          correl <- format(round(correl, 2), nsmall = 2, digits = digits)
          correl[!lower.tri(correl)] <- ""
          print(correl[-1, -p, drop=FALSE], quote = FALSE)
        }
      }
    }
    cat("\n")#- not in S
    invisible(x)
  }

# using qr(<lm>)  as interface to  <lm>$qr :
# #' @title qr.slm
# #'
# #' @param x a numeric or complex matrix whose \code{QR} decomposition is to be computed. Logical matrices are coerced to numeric.
# #' @param ... further arguments passed to or from other methods.
# #'
# #' @seealso
# #' \code{\link[base:qr]{qr}}.
# #'
#' @export
qr.slm <- function(x, ...) {
  if(is.null(r <- x$qr))
    stop("slm object does not have a proper 'qr' component.
         Rank zero or should not have used lm(.., qr=FALSE).")
  r
}

#' @title Predict for slm object
#'
#' @description Predicted values based on \code{slm} object.
#'
#' @param object an object of class \code{slm}.
#' @param newdata an optional data frame in which to look for variables with which to predict.
#'  \code{newdata} must contain only variables and not the intercept.
#'  If omitted, the fitted values are used.
#' @param interval type of interval calculation. It can be only \code{interval = "confidence"}, the default value. It computes
#'  the confidence intervals for \eqn{x' beta}, where \eqn{x'} is a new observation of the design.
#' @param level tolerance/confidence level.
#' @param ... further arguments passed to or from other methods.
#'
#' @details This function produces predicted values, obtained by evaluating the regression function in the frame \code{newdata}
#'  (which defaults to \code{model.frame(object)}). If \code{newdata} is omitted the predictions are based on the data used for the fit.
#'
#' @return This function produces a vector of predictions or a matrix of predictions and bounds with column names \code{fit}, \code{lwr},
#'  and \code{upr} if \code{interval} is set.
#'
#' @importFrom stats predict.lm qnorm lm
#'
#' @seealso
#'  \code{\link[stats:predict.lm]{predict.lm}}.
#'
#' @export
#'
#' @examples
#' data("shan")
#' reg1 = slm(shan$PM_Xuhui ~ . , data = shan, method_cov_st = "fitAR", model_selec = -1)
#' predict(reg1)
#'
#' data("co2")
#' y = as.vector(co2)
#' x = as.vector(time(co2)) - 1958
#' reg2 = slm(y ~ x + I(x^2) + I(x^3) + sin(2*pi*x) + cos(2*pi*x) + sin(4*pi*x) +
#'  cos(4*pi*x) + sin(6*pi*x) + cos(6*pi*x) + sin(8*pi*x) + cos(8*pi*x),
#'  method_cov_st = "fitAR", model_selec = -1)
#' predict(reg2)
predict.slm <- function(object, newdata=NULL, interval="confidence", level=0.95, ...) {
  alpha = 1 - level
  Dn = object@norm_matrix
  Cn = cov_matrix_estimator(object)
  design = cbind(rep(1,length(as.matrix(object$model[1]))),as.matrix(object$model[-1]))
  tdesign = t(design)

  if (interval=="confidence") {
    if (is.null(newdata)) {
      pred = predict.lm(object, interval = interval, level = level, ...)
      #Bhat~N(B,C), x’B~N(xB,x’Cx) IC pour x’B, IC pour une nouvelle observation).
      for (i in seq(1,length(pred[,1]))) {
        pred[i,2] = pred[i,1] - qnorm(1-(alpha/2))*sqrt(design[i,]%*%solve(Dn)%*%Cn%*%solve(Dn)%*%tdesign[,i])
        pred[i,3] = pred[i,1] + qnorm(1-(alpha/2))*sqrt(design[i,]%*%solve(Dn)%*%Cn%*%solve(Dn)%*%tdesign[,i])
      }
      print(pred)
    } else {
      mod = object$model
      names(mod)[1] = "Y"
      reglm = lm(Y~., data=mod)
      pred = predict.lm(reglm, newdata, interval = interval, level = level, ...)
      newdata = cbind(1,as.matrix(newdata))
      tnewdata = t(newdata)
      for (i in seq(1,length(pred[,1]))) {
        pred[i,2] = pred[i,1] - qnorm(1-(alpha/2))*sqrt(newdata[i,]%*%solve(Dn)%*%Cn%*%solve(Dn)%*%tnewdata[,i])
        pred[i,3] = pred[i,1] + qnorm(1-(alpha/2))*sqrt(newdata[i,]%*%solve(Dn)%*%Cn%*%solve(Dn)%*%tnewdata[,i])
      }
      print(pred)
    }
  } else if (interval=="prediction") {
    print("None available")
  } else {
    predict.lm(object, newdata = newdata, ...) #none renvoie la prediction
  }
}


#' @title Plot.slm
#'
#' @description Same function as the \code{\link[stats:plot.lm]{plot.lm}} function.
#'
#' @param x \code{slm} object.
#' @param ... other parameters to be passed through to plotting functions.
#'
#' @return This function returns the graphics of \code{plot.lm(x)}.
#'
#' @importFrom stats lm
#' @importFrom graphics plot
#'
#' @export
#'
#' @examples
#' data("shan")
#' reg = slm(shan$PM_Xuhui ~ . , data = shan, method_cov_st = "fitAR", model_selec = -1)
#' plot(reg)
plot.slm <- function(x, ...) {
  mod = x$model
  names(mod)[1] = "Y"
  reglm = lm(Y~., data=mod)
  plot(reglm, ...)
}

# plot.slm <- function(x, which = c(1:3, 5),
#                      caption = list("Residuals vs Fitted", "Normal Q-Q",
#                                     "Scale-Location", "Cook's distance",
#                                     "Residuals vs Leverage",
#                                     expression("Cook's dist vs Leverage  " * h[ii] / (1 - h[ii]))),
#                      panel = if(add.smooth) panel.smooth else points,
#                      sub.caption = NULL, main = "",
#                      ask = prod(par("mfcol")) < length(which) && dev.interactive(),
#                      ...,
#                      id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
#                      qqline = TRUE, cook.levels = c(0.5, 1.0),
#                      add.smooth = getOption("add.smooth"), label.pos = c(4,2),
#                      cex.caption = 1, cex.oma.main = 1.25) {
#   mod = x$model
#   names(mod)[1] = "Y"
#   reglm = lm(Y~., data=mod)
#
#   plot(reglm, which = c(1:3, 5),
#        caption = list("Residuals vs Fitted", "Normal Q-Q",
#                       "Scale-Location", "Cook's distance",
#                       "Residuals vs Leverage",
#                       expression("Cook's dist vs Leverage  " * h[ii] / (1 - h[ii]))),
#        panel = if(add.smooth) panel.smooth else points,
#        sub.caption = NULL, main = "",
#        ask = prod(par("mfcol")) < length(which) && dev.interactive(),
#        ...,
#        id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
#        qqline = TRUE, cook.levels = c(0.5, 1.0),
#        add.smooth = getOption("add.smooth"), label.pos = c(4,2),
#        cex.caption = 1, cex.oma.main = 1.25)
#
# }

