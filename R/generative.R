#' @title Some stationary processes
#'
#' @description This is a generative function. The user chooses one of the \code{process}: "iid", "AR1", "AR12", "MA12", "Nonmixing", "sysdyn", and it
#'  generates the chosen process. These processes are fully described in the paper of
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @param n sample size.
#' @param process a list of character to choose the process.
#' @param phi a numeric vector with AR parameters if the process is "AR1" or "AR12".
#' @param theta a numeric vector with MA parameters if the process is "MA12".
#'
#' @return This function returns a vector of observations drawn according to the selected process.
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @export
#'
#' @examples
#' generative_process(200,"Nonmixing")
generative_process <- function(n, process = "AR1", phi = "numeric", theta = "numeric"){
  switch(process,
         iid = {out = rt(n,10)^2 - (5/4)},
         # gaussian AR1, phi1=0.7, sigma_2=1
         AR1 = {out = arima.sim(list(ar = c(phi[1])), n)},
         # AR12, phi1=0.5, phi12=0.2, periode = 12, Bruit gaussien sigma_2=1, X_{t} - 0.5 X_{t-1} - 0.2 X_{t-12} = \epsilon_{t}
         AR12 = {out = arima.sim(list(ar=c(phi[1],0,0,0,0,0,0,0,0,0,0,phi[12])), n)},
         # MA12, theta2=0.5, theta3=0.3,theta12=0.2
         MA12 = {out = arima.sim(list(ma=c(0,theta[2],theta[3],0,0,0,0,0,0,0,0,theta[12])), n, rand.gen=rt,df=10)},
         # Nonmixing
         Nonmixing = {out <- array(0,dim=c(1,n))
         out[1,1]=runif(1, min=0, max=1)
         for(i in 2:n) {
           out[1,i]=0.5*(out[1,i-1] + rbinom(1,1,0.5))
         }
         out = qnorm(out, mean=0, sd=5)
         out <- as.vector(out)},
         # dyn. system
         sysdyn = {
           T=function(x,gamma)
           {
             if(x<500000){T<-x+(1/500000)^{gamma}*x^{gamma+1}} else {T<-2*x-1000000}
           }
           gamma = 0.25
           M <- array(0,dim=c(1,n))
           out <- array(0,dim=c(1,n))
           M[1,1] = 51234*runif(1)
           for(i in 2:n){M[1,i]=T(M[1, i-1],gamma)}
           for(i in 1:n){out[1,i]=M[1,i]/1000000}
           out <- as.vector(out)
         }
  )
  return(out)
}

#' @title Some linear model
#'
#' @description This function returns a design for the regression linear model, without the intercept. The user can choose one of the two models:
#'  "mod1" or "mod2". The first model "mod1" contains just one column, equal to \eqn{i^2 + X_i}, \eqn{i=1,...,n}, where \eqn{X} is an AR(1)
#'  process with \code{phi_1 = 0.5}.
#'
#'  The second model "mod2" contains two columns, the first equal to \eqn{log(i) + sin(i) + X_i} and the second equal to \eqn{i}, for \eqn{i=1,...,n}.
#'  The process \eqn{X} is again an AR(1) process with \code{phi_1 = 0.5}. More information about "mod2" is available in the paper of
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm.
#'
#' @param n samples size.
#' @param model a list of character to choose the model.
#'
#' @return This function returns a data-frame which contains a simulated random design.
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#' @export
#'
#' @examples
#' generative_model(500,"mod1")
generative_model <- function(n, model = "mod1"){
  switch(model,
         #1) beta_0 + beta_1(i^2 + X_i), X~AR(1) (cf. article random design, ex.1, modele 1)
         mod1 = {out = list()
         x = arima.sim(list(ar=1/2), n=n, rand.gen=rnorm, sd=3) #AR(1), sd = 3, pour le design
         icarre = seq(1,n)^2
         X = array(1,dim=c(n,1)) #design X
         X[,1] = icarre + x
         X = as.data.frame(X)
         colnames(X) = c("X1")
         out = X
         },
         #2) beta_0 + beta_1(log(i) + sin(i) + X_i) + beta_2 i, X~AR(1) (cf. article random design, ex.1, modele 2)
         mod2 = {out = list()
         x = arima.sim(list(ar=1/2), n=n, rand.gen=rnorm, sd=3) #AR(1), sd = 3, pour le design
         i = seq(1,n)
         X = array(1,dim=c(n,2)) #design X
         X[,1] = log(i) + sin(i) + x
         X[,2] = i
         X = as.data.frame(X)
         colnames(X) = c("X1","X2")
         out = X
         }
  )
  return(out)
}
