QA_SIS <- function(x, y, tau){
  require(quantreg)
  require(splines)
  require(MASS)
  require(survival)
  p = dim(x)[2]
  n = dim(x)[1]
  fit <- numeric(p)
  y <- y - quantile(y, tau)
  x <- scale(x)
  for(j in 1:p){
    x0 <- x[,j]
    knots = quantile(x0, 1:9 / 10)
    a = bs(x0, knots=knots, degree=1)
    b = rq(y~a, tau=tau)
    fit[j] <- sum((b$fitted)^2) / n
  }
  return(fit)
}

## simulated data
fun <- function(X) { return(0.5 * (X[,1] + X[,2] + X[,101]) + matrix(rnorm(dim(X)[1]), dim(X)[1], 1)) }
tau0 = 0.3
n = 500
p = 1000
X = matrix(rnorm(n*p), n, p)
Y <- fun(X)
## Fitted
Fvs <- order(QA_SIS(X, Y, tau0), decreasing=TRUE)
A_num <- max(c(which(Fvs==1), which(Fvs==2), which(Fvs==101)));

