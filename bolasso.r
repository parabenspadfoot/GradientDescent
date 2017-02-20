dev.off
rm(list=ls(all=TRUE))
library("glmnet")

l2norm <- function(x){sum(x*x)^0.5}

############################
# synthetic data generation
############################

p <- 16
r <- 8
n <- 1000
w <- append(rep(0,p-r),rnorm(r, mean = 0, sd = 1))
w <- w/l2norm(w)
w <- runif(1, min = 1/3, max = 1) * w

G <- cbind(replicate(p, rnorm(p)))
Q <- G%*% t(G)
Q <- Q/abs(sum(diag(Q)))

#http://stats.stackexchange.com/questions/120179/generating-data-with-a-given-sample-covariance-matrix
#https://www.r-bloggers.com/simulating-data-following-a-given-covariance-structure/
z <- rnorm(p)
L <-  chol(Q) # indee Q ==  t(L) %*% L
X <-  t(L) %*% matrix(rnorm(p*n), nrow=p, ncol=n)
X <- t(X) # dim(X) is 1000 by 16

sigma <- 0.1 * (sum((X%*%w)^2)/n)^0.5
epsilon <- rnorm(n)*sigma
Y <- X%*%w + epsilon # dim(Y) is 1000 by 1

n <- 800
Xtest <- X[c(801:1000),]
Ytest <- Y[801:1000]
X <- X[c(1:800),]
Y <- Y[c(1:800)]


##############
#  BOLASSO
##############

#max <- 1
#min <- exp(-15)
max <- 1.5*0.0001681229
min <- 0.5*0.0001681229
nLambda <- 10
lambdaVec <- rev(seq(from = min, to = max, by = (max-min)/(nLambda-1))) #as given in graph

#m <- 2 # number of bootstrap replications 
correctness <- c()
m <- c(2,4,8,16,32,64,128,256)
bsWeightsSupportIntersection <-  matrix(rep(TRUE,160),nrow = 16,ncol = 10)
for (i in m){
  bsIndex <- sample(c(1:n),n, replace = TRUE)
  bsX <- X[bsIndex,]
  bsY <- Y[bsIndex]
  bsModel <- glmnet(bsX, bsY, lambda = lambdaVec)
  bsWeights <- bsModel$beta
  bsWeightsSupport <- as.matrix(bsWeights)
  bsWeightsSupport <- as.matrix(bsWeightsSupport > matrix(rep(0.001,160),nrow = 16,ncol = 10)) #arbitrary cutoff at 0.001 by observating bsWeights
  bsWeightsSupportIntersection <- bsWeightsSupport & bsWeightsSupportIntersection 
  
  for (j in c(1:nLambda)){
    capXtest <- Xtest * bsWeightsSupportIntersection[,j]
    capModel <- glmnet(capXtest, Ytest,lambda = lambdaVec[j])
    loadingVec <- capModel$beta
    loadingY <- capXtest %*% loadingVec
    correctness <- append(correctness,sum(sign(Ytest)&sign(loadingY))/200)
  }
}



