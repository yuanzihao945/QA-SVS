DCSISsort=function(X,Y){
  n=dim(X)[1]; ##sample size
  p=dim(X)[2]; ##dimension
  B=matrix(1,n,1);
  C=matrix(1,1,p);
  sxy1=matrix(0,n,p);
  sxy2=matrix(0,n,p);
  sxy3=matrix(0,n,1);
  sxx1=matrix(0,n,p);
  syy1=matrix(0,n,1);
  for (i in 1:n){
    XX1=abs(X-B%*%X[i,]);
    YY1=sqrt(apply((Y-B%*%Y[i])^2,1,sum))
    sxy1[i,]=apply(XX1*(YY1%*%C),2,mean);
    sxy2[i,]=apply(XX1,2,mean);
    sxy3[i,]=mean(YY1);
    XX2=XX1^2;
    sxx1[i,]=apply(XX2,2,mean);
    YY2=YY1^2;
    syy1[i,]=mean(YY2);
  }
  SXY1=apply(sxy1,2,mean);
  SXY2=apply(sxy2,2,mean)*apply(sxy3,2,mean);
  SXY3=apply(sxy2*(sxy3%*%C),2,mean);
  SXX1=apply(sxx1,2,mean);
  SXX2=apply(sxy2,2,mean)^2;
  SXX3=apply(sxy2^2,2,mean);
  SYY1=apply(syy1,2,mean);
  SYY2=apply(sxy3,2,mean)^2;
  SYY3=apply(sxy3^2,2,mean);
  dcovXY=sqrt(SXY1+SXY2-2*SXY3);
  dvarXX=sqrt(SXX1+SXX2-2*SXX3);
  dvarYY=sqrt(SYY1+SYY2-2*SYY3);
  dcorrXY=dcovXY/sqrt(dvarXX*dvarYY);
  A=order(dcorrXY,decreasing=TRUE)
  return (A)
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

p = 5000

for (i in 1:100) {
  
  data <-read.csv(paste('/Users/yuanzihao/Documents/MATLAB/data',
                        i, '.csv', sep = ''), header = FALSE)
  X = as.matrix(data[, 1:p])
  Y1 = matrix(data[, p+1])
  Y2 = matrix(data[, p+2])
  Y3 = matrix(data[, p+3])
  Y4 = matrix(data[, p+4])
  Y5 = matrix(data[, p+5])
  Y6 = matrix(data[, p+6])
  
  Fvs1 <- DCSISsort(X, Y1)
  Fvs2 <- DCSISsort(X, Y2)
  Fvs3 <- DCSISsort(X, Y3)
  Fvs4 <- DCSISsort(X, Y4)
  Fvs5 <- DCSISsort(X, Y5)
  Fvs6 <- DCSISsort(X, Y6)
  
  A_num1[i] = max(c(which(Fvs1==1), which(Fvs1==2), which(Fvs1==101)));
  A_num2[i] = max(c(which(Fvs2==3), which(Fvs2==4), which(Fvs2==102)));
  A_num3[i] = max(c(which(Fvs3==5), which(Fvs3==6), which(Fvs3==103)));
  A_num4[i] = max(c(which(Fvs4==7), which(Fvs4==8), which(Fvs4==104)));
  A_num5[i] = max(c(which(Fvs5==9), which(Fvs5==10), which(Fvs5==105)));
  A_num6[i] = max(c(which(Fvs6==11), which(Fvs6==12), which(Fvs6==106)));
  
}

quantile(A_num1, c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(A_num2, c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(A_num3, c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(A_num4, c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(A_num5, c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(A_num6, c(0.05, 0.25, 0.5, 0.75, 0.95))

write.csv(A_num1, '/Users/yuanzihao/Desktop/DCSIS_2_5000_A_num1.csv')
write.csv(A_num2, '/Users/yuanzihao/Desktop/DCSIS_2_5000_A_num2.csv')
write.csv(A_num3, '/Users/yuanzihao/Desktop/DCSIS_2_5000_A_num3.csv')
write.csv(A_num4, '/Users/yuanzihao/Desktop/DCSIS_2_5000_A_num4.csv')
write.csv(A_num5, '/Users/yuanzihao/Desktop/DCSIS_2_5000_A_num5.csv')
write.csv(A_num6, '/Users/yuanzihao/Desktop/DCSIS_2_5000_A_num6.csv')
