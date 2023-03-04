## Compute Quantile Adaptive sure Independence Screening (QA_SIS)
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


fun1 <- function(X) {
  return(0.5 * (X[,1] + X[,2] + X[,101]) + matrix(rnorm(dim(X)[1]), dim(X)[1], 1))
}


fun2 <- function(X) {
  return(
    0.8 * X[,3] + 0.5 * (X[,4] + 1)^2 + tan(pi * (X[,102] + 1) / 4) +
      matrix(rnorm(dim(X)[1]), dim(X)[1], 1)
  )
}


fun3 <- function(X) {
  return(
    0.5 * exp(X[,5]) + sin(pi * X[,6] / 2) + 
      X[,103] * (X[,103] > quantile(X[,103], 0.6)) +
      matrix(rnorm(dim(X)[1]), dim(X)[1], 1)
  )
}


fun4 <- function(X) {
  return(
    (1 + X[,7] + X[,8])^(-3) * exp(1 + 3 * sin(pi * X[,104] / 2)) +
      matrix(rnorm(dim(X)[1]), dim(X)[1], 1)
  )
}


fun5 <- function(X) {
  return(
    0.5 * X[,9] + tan(pi * (X[,10] - 1) * (X[,105] + 1) / 4) +
      matrix(rnorm(dim(X)[1]), dim(X)[1], 1)
  )
}


fun6 <- function(X) {
  return(
    2 * (X[,11] + 2)^2 * X[,12] * (X[,12] > quantile(X[,12], 0.5)) *
      X[,106] * (X[,106] < quantile(X[,106], 0.5)) +
      matrix(rnorm(dim(X)[1]), dim(X)[1], 1)
  )
}


A_num1 <- matrix(0, 1, 100)
A_num2 <- matrix(0, 1, 100)
A_num3 <- matrix(0, 1, 100)
A_num4 <- matrix(0, 1, 100)
A_num5 <- matrix(0, 1, 100)
A_num6 <- matrix(0, 1, 100)

tau0 = 0.1

for (i in 1:100) {
  X = matrix(rnorm(500*5000), 500, 5000)
  
  Y1 <- fun1(X)
  Y2 <- fun2(X)
  Y3 <- fun3(X)
  Y4 <- fun4(X)
  Y5 <- fun5(X)
  Y6 <- fun6(X)
  
  Fvs1 <- order(QA_SIS(X, Y1, tau0), decreasing=TRUE)
  Fvs2 <- order(QA_SIS(X, Y2, tau0), decreasing=TRUE)
  Fvs3 <- order(QA_SIS(X, Y3, tau0), decreasing=TRUE)
  Fvs4 <- order(QA_SIS(X, Y4, tau0), decreasing=TRUE)
  Fvs5 <- order(QA_SIS(X, Y5, tau0), decreasing=TRUE)
  Fvs6 <- order(QA_SIS(X, Y6, tau0), decreasing=TRUE)
  
  A_num1[i] = max(c(which(Fvs1==1), which(Fvs1==2), which(Fvs1==101)));
  A_num2[i] = max(c(which(Fvs1==3), which(Fvs1==4), which(Fvs1==102)));
  A_num3[i] = max(c(which(Fvs1==5), which(Fvs1==6), which(Fvs1==103)));
  A_num4[i] = max(c(which(Fvs1==7), which(Fvs1==8), which(Fvs1==104)));
  A_num5[i] = max(c(which(Fvs1==9), which(Fvs1==10), which(Fvs1==105)));
  A_num6[i] = max(c(which(Fvs1==11), which(Fvs1==12), which(Fvs1==106)));
  
}

qA_num1 <- quantile(A_num1, c(0.05, 0.25, 0.5, 0.75, 0.95))
qA_num2 <- quantile(A_num2, c(0.05, 0.25, 0.5, 0.75, 0.95))
qA_num3 <- quantile(A_num3, c(0.05, 0.25, 0.5, 0.75, 0.95))
qA_num4 <- quantile(A_num4, c(0.05, 0.25, 0.5, 0.75, 0.95))
qA_num5 <- quantile(A_num5, c(0.05, 0.25, 0.5, 0.75, 0.95))
qA_num6 <- quantile(A_num6, c(0.05, 0.25, 0.5, 0.75, 0.95))
